#include<iomanip>

#include "fit_measurements.h"

using namespace pmp;
using namespace std;


FitMeasure::FitMeasure() {}

FitMeasure::~FitMeasure() {}
//
//
//void FitMeasure::CaculateLaplacianCotMatrix_Test(
//	SurfaceMesh& mesh,
//	Eigen::SparseMatrix<double>& L,
//	std::vector<Tri>& triplets,
//	Eigen::SparseMatrix<double>& b)
//{
//	//std::vector<Tri> triplets;
//	const int p_num = mesh.n_vertices();
//	//std::vector<int> fix{ 275,12478,21,2,3,12222,33,444,555,666,777,6663,4443};
//	//std::vector<int> mov{ 4986 }; //position[-0.316209 0.298498 - 0.088143]
//	std::vector<int> fix{ 1031 };
//	std::vector<int> mov{ 0 };
//	//vertex# 1031
//	//	position[0.997592 0.049009 0.000000]
//	//for (auto v : mesh.vertices())
//	//{
//	//	//圆柱体z坐标小于11固定
//	//	if (mesh.position(v)[2] < -0.8 || mesh.position(v)[2] > 0.6)
//	//	{
//	//		fix.push_back(v.idx());
//	//	}
//	//	//设置z坐标大于18为移动
//	//	//if (mesh.position(v)[2] >= 15)
//	//	//{
//	//	//	mov.push_back(v.idx());
//	//	//}
//	//}
//	int num_fix = fix.size();
//	int num_mov = mov.size();
//	auto fit = mesh.faces_begin();
//	triplets.reserve(20);
//	L.resize(p_num, p_num);
//	do
//	{
//		auto vf = mesh.vertices(*fit);
//		Point p[3];
//		double cot[3];
//		int id[3];
//		for (int i = 0; i < 3; ++i, ++vf)
//		{
//			p[i] = mesh.position(*vf);
//			id[i] = (*vf).idx();
//		}
//		for (int i = 0; i < 3; ++i)
//		{
//			int j = (i + 1) % 3, k = (j + 1) % 3;
//			cot[i] = dot(p[j] - p[i], p[k] - p[i]) / norm(cross(p[j] - p[i], p[k] - p[i]));
//
//			triplets.push_back({ id[j], id[k], -0.5 * cot[i] });
//			triplets.push_back({ id[k], id[j], -0.5 * cot[i] });
//		}
//		for (int i = 0; i < 3; ++i)
//		{
//			triplets.push_back({ id[i], id[i], 0.5*(cot[(i + 1) % 3] + cot[(i + 2) % 3]) });
//		}
//
//	} while (++fit != mesh.faces_end());
//	L.setFromTriplets(triplets.begin(), triplets.end());
//
//	for (int i = 0; i < num_fix; ++i)
//	{
//		triplets.emplace_back(Tri(p_num + i, fix[i], 1));
//	}
//
//	for (int i = 0; i < num_mov; ++i)
//	{
//		triplets.emplace_back(Tri(p_num + num_fix + i, mov[i], 1));
//	}
//	Eigen::SparseMatrix<double> v;
//	v.resize(p_num, 3);
//	std::vector<Tri> tri_v;
//	tri_v.reserve(3);
//	int i = 0;
//	for (const auto &v_it : mesh.vertices())
//	{
//		Point t = mesh.position(v_it);
//		tri_v.emplace_back(Tri(i, 0, t[0]));
//		tri_v.emplace_back(Tri(i, 1, t[1]));
//		tri_v.emplace_back(Tri(i, 2, t[2]));
//		i++;
//	}
//	v.setFromTriplets(tri_v.begin(), tri_v.end());
//	Eigen::SparseMatrix<double> A;
//	A.resize(p_num + num_fix + num_mov, p_num);
//	A.setFromTriplets(triplets.begin(), triplets.end());
//	b = L * v;
//	b.conservativeResize(p_num + num_fix + num_mov, 3);
//	for (int i = 0; i < num_fix; ++i)
//	{
//		for (int j = 0; j < 3; ++j)
//		{
//			b.insert(p_num + i, j) = v.coeff(fix[i], j);
//		}
//	}
//
//	//b(3 * p_num + 3 * num_fix + 0) = -0.55;
//	//b(3 * p_num + 3 * num_fix + 1) = 0.35;
//	//b(3 * p_num + 3 * num_fix + 2) = -0.088;
//	//for (int i = 0; i < num_mov; ++i)
//	//{
//	//	for (int j = 0; j < 3; ++j)
//	//	{
//	//		if (j == 2)
//	//			b.coeffRef(p_num + num_fix + i, j) = v.coeff(mov[i], j) + 5.0;
//	//		else
//	//			b.coeffRef(p_num + num_fix + i, j) = v.coeff(mov[i], j);
//	//	}
//	//}
//	b.insert(p_num + num_fix, 0) = v.coeff(mov[0], 0) - 15;
//	b.insert(p_num + num_fix, 1) = v.coeff(mov[0], 1);
//	b.insert(p_num + num_fix, 2) = v.coeff(mov[0], 2);
//	//std::cout << b << std::endl;
//	//position[-0.316209 0.298498 - 0.088143]
//
//	auto AT = A.transpose();
//	ShowMessage(string("AT"));
//
//	Eigen::SparseMatrix<double> ATA = AT * A;
//	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(ATA);
//	ShowMessage(string("ATA"));
//
//	if (solver.info() != Eigen::Success)
//		ShowMessage(string(">Solve Failed"));
//	//Eigen::Matrix3Xi facets;
//	//binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "facets").c_str(), facets);
//	Eigen::MatrixXd new_vertice(p_num, 3);
//	new_vertice = solver.solve(AT*b);
//	//new_vertice.col(1) = solver.solve(AT*b.col(1));
//	//new_vertice.col(2) = solver.solve(AT*b.col(2));
//	ShowMessage(string(">Solve Success"));
//
//	//meshio::SaveObj((BIN_DATA_PATH + "res.obj").c_str(), res_verts, facets);
//
//	for (auto v : mesh.vertices())
//	{
//		Point& pos = mesh.position(v);
//		pos[0] = new_vertice.coeff(v.idx(), 0);
//		pos[1] = new_vertice.coeff(v.idx(), 1);
//		pos[2] = new_vertice.coeff(v.idx(), 2);
//	}
//	mesh.write((BIN_DATA_PATH + "res.obj").c_str());
//}
//
//
//void FitMeasure::ConstructCoefficientMatrixBottom_Test(
//	SurfaceMesh& mesh,
//	std::vector<std::vector<int>>& point_idx,
//	const Eigen::MatrixXd& measurements,
//	Eigen::SparseMatrix<double>& b,
//	const Eigen::MatrixXd& input_m,
//	std::vector<Tri>& triplets)
//{
//	vector<Point> Vertice;
//	for (const auto &v_it : mesh.vertices())
//	{
//		Point t = mesh.position(v_it);
//		Vertice.push_back(t);
//	}
//	num_edge_all = 3;
//	vector<int> point_idx_t{ 12478, 275 };
//	vector<int> fix{ 10941, 11120 };
//	//std::vector<Tri> triplets;
//	num_verts = mesh.n_vertices();
//	b.resize(3, 3);
//	triplets.reserve(7);
//	int row = 0;
//	//遍历每个尺寸的每条边
//	for (int j = 0; j < 1; ++j)
//	{
//		int edge_0 = point_idx_t[0];
//		int edge_1 = point_idx_t[1];
//		const Eigen::Vector3d v0 = Vertice[edge_0];
//		const Eigen::Vector3d v1 = Vertice[edge_1];
//		const Eigen::Vector3d edge_01 = v1 - v0;
//		const double edge_len = edge_01.norm();
//		double target_len = 1.8;
//		//std::cout << edge_len << std::endl;
//		Eigen::Vector3d auxd = (edge_01 / edge_len) * target_len;
//
//		int _row = num_verts + j;
//		int _col0 = edge_0;
//		int _col1 = edge_1;
//		triplets.emplace_back(Tri(_row, _col0, -1));
//		triplets.emplace_back(Tri(_row, _col1, 1));
//		for (int k = 0; k < 3; ++k)
//		{
//			b.insert(j, k) = auxd[k];
//		}
//		//std::cout << "row: " << row + j * 3 + k << " " << "col: " << edge_0 * 3 + k << " " << auxd[k] << std::endl;
//	}
//
//	for (int i = 0; i < fix.size(); ++i)
//	{
//		triplets.emplace_back(Tri(num_verts + i, fix[i], 1));
//	}
//
//	for (int i = 0; i < fix.size(); ++i)
//	{
//		for (int j = 0; j < 3; ++j)
//		{
//
//			//b.insert(1 + i, j) = vertices.coeff(j, fix[i]);
//		}
//	}
//}
//
//void FitMeasure::FitMeasurements_Test(
//	SurfaceMesh& mesh,
//	Eigen::Matrix3Xd& res_verts,
//	const Eigen::SparseMatrix<double>& L,
//	const Eigen::Matrix3Xd& vertices,
//	Eigen::SparseMatrix<double> & b2,
//	std::vector<std::vector<int>>& point_idx,
//	Eigen::SparseMatrix<double> A)
//{
//	Eigen::SparseMatrix<double> v(num_verts, 3);
//	Mat2Vec(v, vertices);
//
//	Eigen::SparseMatrix<double> b1 = L * v;
//	ShowMessage(string("b1 = L * v"));
//
//	Eigen::SparseMatrix<double> b(num_verts + num_edge_all, 3);
//	//ConstructB(b, b1, b2);
//
//	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
//	auto AT = A.transpose();
//	ShowMessage(string("AT"));
//
//	Eigen::SparseMatrix<double> ATA = AT * A;
//	solver.compute(ATA);
//	ShowMessage(string("ATA"));
//
//	if (solver.info() != Eigen::Success)
//		ShowMessage(string(">Solve Failed"));
//
//	Eigen::SparseMatrix<double> new_vertice = solver.solve(AT * b);
//	ShowMessage(string(">Solve Success"));
//	//res_verts = Eigen::Map<Eigen::Matrix3Xd>(vecV.data(), vertices.cols(), 3);
//	for (auto v : mesh.vertices())
//	{
//		Point& pos = mesh.position(v);
//		pos[0] = new_vertice.coeff(v.idx(), 0);
//		pos[1] = new_vertice.coeff(v.idx(), 1);
//		pos[2] = new_vertice.coeff(v.idx(), 2);
//	}
//	mesh.write((BIN_DATA_PATH + "res.obj").c_str());
//}
//void FitMeasure::Mat2Vec_Test(
//	Eigen::SparseMatrix<double>& v,
//	const std::vector<Point>& Vertice)
//{
//	for (int i = 0; i < num_verts; ++i)
//	{
//		for (int j = 0; j < 3; ++j)
//		{
//			v.coeffRef(i, j) = Vertice[i][0];
//		}
//	}
//}
//

//void CalcEnergy(double& enyg, Eigen::Matrix3Xd& vertices)
//{
//	for (size_t i_ = 1, i = 1; i_ < point_idx.size(); ++i_)
//	{
//		//剔除三个围长
//		if (i_ == 2 || i_ == 3 || i_ == 4)
//		{
//			continue;
//		}
//		size_t n = point_idx[i_].size();
//		for (size_t j = 0; j < n - 1; ++j)
//		{
//			int id1 = point_idx[i_][j], id2 = point_idx[i_][j + 1];
//			pmp::vec3 p1; p1[0] = vertices.coeff(0, id1); p1[1] = vertices.coeff(1, id1); p1[2] = vertices.coeff(2, id1);
//			pmp::vec3 p2; p2[0] = vertices.coeff(0, id2); p2[1] = vertices.coeff(1, id2); p2[2] = vertices.coeff(2, id2);
//			double cur_len = distance(p1, p2);
//			double p = std::pow(cur_len, 2);
//			double l = std::pow(CalcTargetLen(measurements, cur_len, i), 2);
//			enyg = enyg + std::pow(p - l, 2);
//			std::cout << setprecision(15);
//			//std::cout << std::fixed << p << "  " << l << "  " << enyg << std::endl;
//		}
//		i++;
//	}
//}

/*!
*@brief  计算欧式距离的梯度
*@param[out]
*@param[in]  Eigen::MatrixXd & gradient
*@param[in]  SurfaceMesh & mesh
*@param[in]  std::vector<int> point_idx
*@param[in]  Eigen::MatrixXd & input_measure
*@return     void
*/
//void FitMeasure::CalcEuclideanGradient(Eigen::VectorXd & gradient, Matrix3Xd & vertices)
//{
//	pmp::vec3 grad;
//	int id1 = point_idx[0][0], id2 = point_idx[0][1];
//	pmp::vec3 p1; p1[0] = vertices.coeff(0, id1); p1[1] = vertices.coeff(1, id1); p1[2] = vertices.coeff(2, id1);
//	pmp::vec3 p2; p2[0] = vertices.coeff(0, id2); p2[1] = vertices.coeff(1, id2); p2[2] = vertices.coeff(2, id2);
//	grad = 4 * (std::pow(distance(p1, p2), 2) - std::pow(input_m.coeff(0, 0), 2))*(p1 - p2);
//	for (int i = 0; i < 3; ++i)
//	{
//		gradient(id1 * 3 + i) = grad[i];
//		gradient(id2 * 3 + i) = -grad[i];
//	}
//}

/*!
*@brief  计算测地距离的梯度
*@param[out]
*@param[in]  Eigen::VectorXd & gradient
*@param[in]  Matrix3Xd vertices
*@param[in]  Eigen::MatrixXd measurements
*@return     void
*/
//void FitMeasure::CalcGeodesicGradient(Eigen::VectorXd & gradient, Matrix3Xd vertices, Eigen::MatrixXd & measurements, Eigen::MatrixXd & input_m)
//{
//	pmp::vec3 grad;
//	//从1开始因为poin_idx[0]是欧式距离
//	for (size_t i_ = 1, i = 1; i_ < point_idx.size(); ++i_)
//	{
//		if (i_ == 2 || i_ == 3 || i_ == 4)
//		{
//			continue;
//		}
//		size_t n = point_idx[i_].size();
//		for (size_t j = 0; j < n - 1; ++j)
//		{
//			int id1 = point_idx[i_][j], id2 = point_idx[i_][j + 1];
//			//std::cout << id1 << "  " << id2 << std::endl;
//			pmp::vec3 p1; p1[0] = vertices.coeff(0, id1); p1[1] = vertices.coeff(1, id1); p1[2] = vertices.coeff(2, id1);
//			pmp::vec3 p2; p2[0] = vertices.coeff(0, id2); p2[1] = vertices.coeff(1, id2); p2[2] = vertices.coeff(2, id2);
//			float cur_len = distance(p1, p2);
//			grad = 4 * (std::pow(cur_len, 2) - std::pow(CalcTargetLen(measurements, cur_len, i, input_m), 2))*(p1 - p2);
//			for (size_t k = 0; k < 3; ++k)
//			{
//				gradient(id1 * 3 + k) += grad[k];
//				gradient(id2 * 3 + k) += -grad[k];
//			}
//		}
//		++i;
//	}
//}


/*!
*@brief  计算目标边长
*@param[out]
*@param[in]  Eigen::MatrixXd measurements
*@param[in]  const float cur_len  当前网格上这条边长
*@param[in]  int index  第index个尺寸
*@return     float
*/
float FitMeasure::CalcTargetLen(
	const Eigen::MatrixXd& one_measure,
	const float& cur_len,
	const int& index,
	const Eigen::MatrixXd& input_m)
{
	float m = one_measure.coeff(index + 1, 0) / 1000.0;
	float cur_len_p = cur_len / m;
	float target_len_p = input_m.coeff(index, 0) / 1000.0;
	float target_len = cur_len_p * target_len_p;
	return target_len;
}

/*!
*@brief  计算拉普拉斯权值矩阵
*@param[out] 拉普拉斯稀疏矩阵
*@param[in]  const SurfaceMesh & mesh  待求拉普拉斯矩阵的原始网格
*@param[in]  Eigen::SparseMatrix<float> & L  拉普拉斯矩阵，是个大型稀疏矩阵
*@return     void
*/
void FitMeasure::CaculateLaplacianCotMatrix(
	const SurfaceMesh& mesh)
{
	//std::vector<Tri> triplets;
	const int p_num = mesh.n_vertices();
	L.resize(p_num, p_num);
	auto fit = mesh.faces_begin();
	do {
		auto vf = mesh.vertices(*fit);
		Point p[3];
		int id[3];
		for (int i = 0; i < 3; ++i, ++vf)
		{
			p[i] = mesh.position(*vf);
			id[i] = (*vf).idx();
		}
		SetTriplets(p, id);
	} while (++fit != mesh.faces_end());
	L.setFromTriplets(triplets_A.begin(), triplets_A.end());
	std::cout << "Calc Laplacian Done !" << std::endl;
}

void FitMeasure::SetTriplets(
	vec3 p[3],
	int id[3])
{
	triplets_A.reserve(7);
	float cot[3];
	for (int i = 0; i < 3; ++i)
	{
		int j = (i + 1) % 3, k = (j + 1) % 3;
		cot[i] = dot(p[j] - p[i], p[k] - p[i]) / norm(cross(p[j] - p[i], p[k] - p[i]));
		triplets_A.push_back({ id[j],id[k], -0.5 * cot[i] });
		triplets_A.push_back({ id[k],id[j], -0.5 * cot[i] });

	}
	for (int i = 0; i < 3; ++i)
	{
		triplets_A.push_back({ id[i], id[i], 0.5*(cot[(i + 1) % 3] + cot[(i + 2) % 3]) });
	}

}

/*!
*@brief  保存每个尺寸边的数量
*@param[out]
*@param[in]  std::vector<std::vector<int>> & point_idx
*@return     void
*/
void FitMeasure::SaveEdge(
	std::vector<std::vector<int>>& point_idx)
{
	//保存每个尺寸的边数量
	//std::vector<int> edge;//size = 18
	for (auto& num_v : point_idx)
	{
		if (num_v.size() > 2)
		{
			edge.push_back(num_v.size());
		}
		else
		{
			edge.push_back(1);
			//continue;
		}
	}
}

void FitMeasure::ConstructCoefficientMatrixBottom()
{
	C.resize(num_edge_all, num_verts);
	C.setFromTriplets(triplets_C.begin(), triplets_C.end());
}
/*!
*@brief  构建系数矩阵，A矩阵下半部分
*@param[out]
*@param[in]  Eigen::SparseMatrix<double> & A
*@param[in]  std::vector<std::vector<int>> & point_idx
*@param[in]  const Eigen::Matrix3Xd & vertices
*@param[in]  const Eigen::MatrixXd & measurements
*@param[in]  Eigen::VectorXd & b
*@param[in]  Eigen::MatrixXd & input_m
*@return     void
*/

void FitMeasure::ConstructCoefficientMatrix(
	std::vector<std::vector<int>>& point_idx,
	const Eigen::Matrix3Xd &vertices,
	const Eigen::MatrixXd& one_measure,
	const Eigen::MatrixXd& input_m)
{
	num_verts = vertices.cols();
	num_measure = point_idx.size();	//18
	SaveEdge(point_idx);
	num_edge_all = std::accumulate(edge.begin(), edge.end(), 0);	//所有边的数量1246
	b_down.resize(num_edge_all, 3);
	SetTriplets(vertices, input_m, point_idx, one_measure);
	A.resize(num_verts + num_edge_all, num_verts);
	A.setFromTriplets(triplets_A.begin(), triplets_A.end());
	ShowMessage(string("Construct A"));
}

void FitMeasure::SetTriplets(
	const Eigen::Matrix3Xd& vertices,
	const Eigen::MatrixXd& input_m,
	const std::vector<std::vector<int>>& point_idx,
	const Eigen::MatrixXd& one_measure)
{
	triplets_A.reserve(7);
	int row = 0;
	for (int i = 0; i < num_measure; ++i)
	{
		//遍历每个尺寸的每条边
		for (int j = 0; j < edge[i]; ++j)
		{
			int edge_0;
			int edge_1;
			//i == 0身高信息
			if (i == 0)
			{
				edge_0 = point_idx[i][j];
				edge_1 = point_idx[i][j + 1];
			}
			else
			{
				edge_0 = point_idx[i][j % edge[i]];
				edge_1 = point_idx[i][(j + 1) % edge[i]];
			}
			const Eigen::Vector3d v0 = vertices.col(edge_0);
			const Eigen::Vector3d v1 = vertices.col(edge_1);
			const Eigen::Vector3d edge_01 = v1 - v0;
			const double edge_len = edge_01.norm();
			double target_len = CalcTargetLen(one_measure, edge_len, i, input_m);
			//std::cout << edge_len << std::endl;
			Eigen::Vector3d auxd = (edge_01 / edge_len) * target_len;

			int _row = num_verts + row + j;
			int _col0 = edge_0;
			int _col1 = edge_1;
			triplets_A.emplace_back(Tri(_row, _col0, -1));
			triplets_A.emplace_back(Tri(_row, _col1, 1));
			triplets_C.emplace_back(Tri(row + j, _col0, -1));
			triplets_C.emplace_back(Tri(row + j, _col1, 1));
			for (int k = 0; k < 3; ++k)
			{
				b_down.insert(row + j, k) = auxd[k];
			}
			//std::cout << "row: " << row + j * 3 + k << " " << "col: " << edge_0 * 3 + k << " " << auxd[k] << std::endl;
		}
		row += edge[i];
	}
	//std::cout << row << std::endl;
}

/*!
*@brief  	将顶点矩阵大小3*V转成3V*1
*@param[out]
*@param[in]  Eigen::SparseVector<double> & v
*@param[in]  const Eigen::Matrix3Xd & vertices
*@return     void
*/
void FitMeasure::Mat2Vec(
	Eigen::SparseMatrix<double>& v,
	const Eigen::Matrix3Xd& vertices)
{
	for (int i = 0; i < num_verts; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			v.coeffRef(i, j) = vertices.coeff(j, i);
		}
	}
}

/*!
*@brief  构建矩阵b，即方程右边
*@param[out]
*@param[in]  Eigen::SparseMatrix<double> & b1
*@param[in]  Eigen::VectorXd & b2
*@return     void
*/
void FitMeasure::ConstructB(
	Eigen::SparseMatrix<double>& b1,
	Eigen::SparseMatrix<double>& b2)
{
	b.resize(num_verts + num_edge_all, 3);
	for (int i = 0; i < num_verts; ++i)
		for (int j = 0; j < 3; ++j)
		{
			b.insert(i, j) = b1.coeff(i, j);
		}

	for (int i = 0; i < num_edge_all; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			b.insert(num_verts + i, j) = b2.coeff(i, j);
		}
	}
	ShowMessage(string("FillB"));
}


void FitMeasure::SaveObj(
	SurfaceMesh& mesh,
	Eigen::SparseMatrix<double>& new_vertice)
{
	for (auto v : mesh.vertices())
	{
		Point& pos = mesh.position(v);
		pos[0] = new_vertice.coeff(v.idx(), 0);
		pos[1] = new_vertice.coeff(v.idx(), 1);
		pos[2] = new_vertice.coeff(v.idx(), 2);
	}
	mesh.write((BIN_DATA_PATH + "res.obj").c_str());
}
void FitMeasure::SaveObj(
	SurfaceMesh& mesh,
	Eigen::MatrixX3d& new_vertice)
{
	for (auto v : mesh.vertices())
	{
		Point& pos = mesh.position(v);
		pos[0] = new_vertice.coeff(v.idx(), 0);
		pos[1] = new_vertice.coeff(v.idx(), 1);
		pos[2] = new_vertice.coeff(v.idx(), 2);
	}
	mesh.write((BIN_DATA_PATH + "res.obj").c_str());
}

void FitMeasure::ShowMessage(const string& msg)
{
	std::cout << msg.c_str() << std::endl;
}



//将归一化的尺寸恢复为原来大小
void FitMeasure::RecoverMeasure(Eigen::MatrixXd& measurelist, Eigen::MatrixXd& one_measure)
{
	Eigen::MatrixXd mean_measure, std_measure;
	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "mean_measure").c_str(), mean_measure);
	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "std_measure").c_str(), std_measure);
	one_measure = one_measure.cwiseProduct(std_measure);
	one_measure += mean_measure;
}

Eigen::MatrixX3d FitMeasure::AULSolver(
	const Eigen::Matrix3Xd& vertices)
{
	Eigen::MatrixX3d res_vert; res_vert.resize(num_verts, 3);
	for (int i = 0; i < 3; ++i)
	{
		Eigen::MatrixXd vertices_one = vertices.row(i);
		b_col = i;
		res_vert.col(i) = SolveOneDimension(vertices_one);
	}
	return res_vert;
}

Eigen::MatrixXd FitMeasure::SolveOneDimension(
	Eigen::MatrixXd& vertices_one)
{
	vertices_one.resize(num_verts, 1);
	auto CT = C.transpose();
	Eigen::MatrixXd vertices_one_next = vertices_one;
	Eigen::MatrixXd alpha_next(alpha);
	alpha.resize(num_edge_all, 1); alpha.setConstant(0.1);

	do
	{
		alpha = alpha_next;
		vertices_one_next = LBFGS(vertices_one);
		alpha_next = alpha + beta * (C * vertices_one_next - b_down.col(b_col));
	} while (vertices_one.norm() - vertices_one_next.norm() < eps);
	return vertices_one_next;
}

Eigen::SparseMatrix<double> FitMeasure::CalcGradient(
	Eigen::MatrixXd& V0)
{
	Eigen::MatrixXd F_up = L * V0;
	Eigen::MatrixXd F_down = C * V0;
	auto LT = L.transpose();
	auto CT = C.transpose();
	Eigen::SparseMatrix<double> grad_fx = 2 * LT * (L * V0 - b_up.col(b_col) + b_up.col(b_col));
	return grad_fx;
}

void FitMeasure::SetGrad(Eigen::MatrixXd& grad_t, real_1d_array& grad)
{
	for (int i = 0; i < grad_t.size(); ++i)
	{
		grad[i] = grad_t.coeff(i, 0);
	}
}

Eigen::MatrixXd FitMeasure::LBFGS(Eigen::MatrixXd& vertices_one)
{
	alglib::real_1d_array x;
	x.attach_to_ptr(num_verts, vertices_one.data());
	Eigen::VectorXd scale;
	scale.resize(num_verts);
	scale.setOnes();
	real_1d_array s;
	s.setcontent(num_verts, scale.data());

	double epsg = 0;
	double epsf = 0;
	double epsx = 1e-8;
	ae_int_t maxits = 0;
	minlbfgsstate state;
	minlbfgscreate(1, x, state);
	minlbfgssetcond(state, epsg, epsf, epsx, maxits);
	minlbfgssetscale(state, s);


	minlbfgsoptguardsmoothness(state);
	minlbfgsoptguardgradient(state, 0.001);
	minlbfgsreport rep;

	alglib::minlbfgsoptimize(state, function1_grad);
	minlbfgsresults(state, x, rep);

	Eigen::MatrixXd res;
	array2mat(x, res);

	return res;

	//
	// Check that OptGuard did not report errors
	//
	// NOTE: want to test OptGuard? Try breaking the gradient - say, add
	//       1.0 to some of its components.
	//
	optguardreport ogrep;
	minlbfgsoptguardresults(state, ogrep);
	printf("%s\n", ogrep.badgradsuspected ? "true" : "false"); // EXPECTED: false
	printf("%s\n", ogrep.nonc0suspected ? "true" : "false"); // EXPECTED: false
	printf("%s\n", ogrep.nonc1suspected ? "true" : "false"); // EXPECTED: false
}

void FitMeasure::array2mat(real_1d_array& x, Eigen::MatrixXd& res)
{
	res.resize(num_verts, 1);
	for (int i = 0; i < num_verts; ++i)
	{
		res.coeffRef(i, 0) = x[i];
	}
}

void function1_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr)
{
	FitMeasure fit;
	int num_verts = fit.num_verts;
	int num_edge_all = fit.num_edge_all;
	alglib::real_1d_array temp(x);
	Eigen::MatrixXd V0 = Eigen::Map<Eigen::MatrixXd>(temp.getcontent(), num_verts, 1);
	Eigen::MatrixXd F_up = fit.L * V0;
	Eigen::MatrixXd F_down = fit.C * V0;
	auto CT = fit.C.transpose();
	Eigen::MatrixXd alpha = fit.alpha;
	alpha.resize(num_edge_all, 1); alpha.setConstant(0.1);
	double beta = fit.beta;

	for (int i = 0; i < num_verts; ++i)
	{
		auto t = F_up.coeff(i, 0) - fit.b_up.coeff(i, fit.b_col);
		func += t;
	}
	for (int i = 1; i <= num_edge_all; ++i)
	{
		auto t = F_down.coeff(i - 1, 0) - fit.b_down.coeff(i - 1, fit.b_col);
		func += t;
	}

	Eigen::MatrixXd grad_t = fit.CalcGradient(V0) + CT * alpha +
		beta * CT*(fit.C * V0 - fit.b_down.col(fit.b_col));
	fit.SetGrad(grad_t, grad);
}


void FitMeasure::ErrorAnalysis(
	Measure	measure,
	std::vector<std::vector<std::vector<double>>>& control_points,
	Eigen::MatrixXd input_m)
{
	Eigen::Matrix3Xd vertex_1;
	Eigen::Matrix3Xd vertex_2;
	Eigen::Matrix3Xi facets;
	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "facets").c_str(), facets);
	meshio::ReadObj((DATASET_PATH + "R5.obj").c_str(), vertex_1, facets);
	meshio::ReadObj((BIN_DATA_PATH + "res.obj").c_str(), vertex_2, facets);

	Eigen::MatrixXd measure_1;
	Eigen::MatrixXd measure_2;
	Eigen::MatrixXd error;
	measure_1.resize(19, 1);
	measure_2.resize(19, 1);
	measure_1 = measure.CalcMeasure(control_points, vertex_1, facets);
	measure_2 = measure.CalcMeasure(control_points, vertex_2, facets);
	std::cout << setprecision(15);
	error = input_m - measure_2.middleRows(1, 18);
	cout << std::fixed << input_m << endl;
	cout << "-----------" << endl;
	cout << std::fixed << measure_2 << endl;
	cout << "-----------" << endl;
	cout << std::fixed << error << endl;
}

/*!
*@brief  拟合尺寸
*@param[out]
*@param[in]  Eigen::Matrix3Xd & res_verts
*@param[in]  Eigen::SparseMatrix<double> & C
*@param[in]  const Eigen::Matrix3Xd & vertices
*@param[in]  Eigen::VectorXd & b2
*@param[in]  std::vector<std::vector<int>> & point_idx
*@return     void
*/
void FitMeasure::FitMeasurements(
	SurfaceMesh& mesh,
	const Eigen::Matrix3Xd& vertices,
	const Eigen::MatrixXd& one_measure,
	std::vector<std::vector<int>>& point_idx,
	const Eigen::MatrixXd& input_m)
{
	CaculateLaplacianCotMatrix(mesh);
	ConstructCoefficientMatrix(point_idx, vertices, one_measure, input_m);
	ConstructCoefficientMatrixBottom();

	Eigen::SparseMatrix<double> v(num_verts, 3);
	Mat2Vec(v, vertices);

	Eigen::SparseMatrix<double> b1 = L * v;
	b_up = b1;
	ShowMessage(string("b1 = L * v"));
	ConstructB(b1, b_down);
	//Eigen::MatrixX3d new_vertice = AULSolver(vertices);
	//SaveObj(mesh, new_vertice);
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	auto AT = A.transpose();
	ShowMessage(string("AT"));

	Eigen::SparseMatrix<double> ATA = AT * A;
	solver.compute(ATA);
	ShowMessage(string("ATA"));

	if (solver.info() != Eigen::Success)
		ShowMessage(string(">Solve Failed"));

	Eigen::SparseMatrix<double> new_vertice = solver.solve(AT * b);
	ShowMessage(string(">Solve Success"));
	SaveObj(mesh, new_vertice);
}
