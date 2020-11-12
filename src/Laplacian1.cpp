
#include "Laplacian1.h"

using namespace std;
using namespace Eigen;
using namespace pmp;

LaplaceDeformation::LaplaceDeformation()
{
}

//初始化：设置固定锚点和移动锚点已经移动锚点移动后的坐标
void LaplaceDeformation::InitializeMesh()
{
	move_anchor_coord.resize(1);
	move_anchor_coord[0][0] = -0.245001;
	move_anchor_coord[0][1] = 0.8;
	move_anchor_coord[0][2] = -0.062755;
}

//主函数
void LaplaceDeformation::AllpyLaplaceDeformation()
{
	//获取数据
	SurfaceMesh mesh;
	mesh.read((DATASET_PATH + "AVE.obj").c_str());
	cout << "获取数据完成" << endl;

	if (mesh.n_vertices() == 0)
		return;
	InitializeMesh();
	cout << "初始化完成" << endl;
	std::vector<T> triplets;
	Compute_CotMatrix(mesh, triplets);
	cout << "余切矩阵计算完成！" << endl;

	BuildATtimesAMatrix(mesh, triplets);
	cout << "ATA" << endl;

	BuildATtimesbMatrix(mesh);
	cout << "ATb" << endl;

	//X = ATA.ldlt().solve(B) //X: n*3, B:N*3
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(ATA);
	new_vertice = chol.solve(ATb);
	cout << "求解完成" << endl;

	SetNewcord(mesh);
	cout << "更新" << endl;

	mesh.write((BIN_DATA_PATH + "res.obj").c_str());
}


void LaplaceDeformation::Compute_CotMatrix(const SurfaceMesh& mesh, std::vector<T>& triplets)
{
	auto pts = mesh.get_vertex_property<Point>("v:point");
	auto fit = mesh.faces_begin();
	triplets.reserve(20);
	const int p_num = mesh.n_vertices();
	L.resize(p_num, p_num);

	do
	{
		auto vf = mesh.vertices(*fit);
		Point p[3];
		double cot[3];
		int id[3];
		for (int i = 0; i < 3; ++i, ++vf)
		{
			p[i] = pts[*vf];
			id[i] = (*vf).idx();
		}

		for (int i = 0; i < 3; ++i)
		{
			int j = (i + 1) % 3, k = (j + 1) % 3;
			cot[i] = dot(p[j] - p[i], p[k] - p[i]) / norm(cross(p[j] - p[i], p[k] - p[i]));

			triplets.push_back({ id[j], id[k], -0.5 * cot[i] });
			triplets.push_back({ id[k], id[j], -0.5 * cot[i] });
		}
		for (int i = 0; i < 3; ++i)
		{
			triplets.push_back({ id[i], id[i], 0.5*(cot[(i + 1) % 3] + cot[(i + 2) % 3]) });
		}

	} while (++fit != mesh.faces_end());
	L.setFromTriplets(triplets.begin(), triplets.end());

}

void LaplaceDeformation::BuildATtimesAMatrix(const SurfaceMesh & mesh, std::vector<T>& triplets)
{
	int n_fix_anchors = fix_anchor_idx.size(), n_move_anchors = move_anchor_idx.size(), points_num = mesh.vertices_size();
	for (int i = 0; i < n_fix_anchors; ++i)
	{
		triplets.emplace_back(T(points_num + i, fix_anchor_idx[i], 1));
	}

	for (int i = 0; i < n_move_anchors; ++i)
	{
		triplets.emplace_back(T(points_num + n_fix_anchors + i, move_anchor_idx[i], 1));
	}
	A.resize(points_num + n_fix_anchors + n_move_anchors, points_num);
	A.setFromTriplets(triplets.begin(), triplets.end());
	//printf("Laplace 矩阵计算完成\n");
	//cout << A << endl;
	AT = A.transpose();
	//cout << AT << endl;
	ATA = AT * A;
	//cout << ATA << endl;
}

void LaplaceDeformation::BuildATtimesbMatrix(const SurfaceMesh & mesh)
{
	int n_fix_anchors = fix_anchor_idx.size(), n_move_anchors = move_anchor_idx.size(), points_num = mesh.vertices_size();
	SparseMatrix<double> v(points_num, 3);
	std::vector<T> tripletlist;
	tripletlist.reserve(3);
	int i = 0;
	vector<Point> Vertice;

	for (const auto &v_it : mesh.vertices())
	{
		Point t = mesh.position(v_it);
		Vertice.push_back(t);
		tripletlist.push_back(T(i, 0, t[0]));
		tripletlist.push_back(T(i, 1, t[1]));
		tripletlist.push_back(T(i, 2, t[2]));
		i++;
	}
	v.setFromTriplets(tripletlist.begin(), tripletlist.end());

	// 根据laplace矩阵计算出所有点的的laplace坐标
	b = L * v;
	b.conservativeResize(points_num + n_fix_anchors + n_move_anchors, 3);

	// 用形变前坐标对固定锚点坐标进行赋值
	for (auto i = 0; i < n_fix_anchors; i++)
	{

		b.insert(i + points_num, 0) = Vertice[fix_anchor_idx[i]][0];
		b.insert(i + points_num, 1) = Vertice[fix_anchor_idx[i]][1];
		b.insert(i + points_num, 2) = Vertice[fix_anchor_idx[i]][2];
	}

	// 用形变后坐标对移动锚点坐标进行赋值
	for (auto i = 0; i < n_move_anchors; i++)
	{
		b.insert(i + points_num + n_fix_anchors, 0) = move_anchor_coord[i][0];
		b.insert(i + points_num + n_fix_anchors, 1) = move_anchor_coord[i][1];
		b.insert(i + points_num + n_fix_anchors, 2) = move_anchor_coord[i][2];
	}
	// 计算三个轴上的 ATb 向量
	ATb = AT * b;
}

//更新坐标
void LaplaceDeformation::SetNewcord(SurfaceMesh & mesh)
{
	for (auto v : mesh.vertices())
	{
		Point& pos = mesh.position(v);
		pos[0] = new_vertice.coeff(v.idx(), 0);
		pos[1] = new_vertice.coeff(v.idx(), 1);
		pos[2] = new_vertice.coeff(v.idx(), 2);
	}
}
