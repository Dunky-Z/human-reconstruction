#include<iomanip>

#include "fit_measurements.h"

using namespace pmp;
using namespace std;




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
float FitMeasure::CalcTargetLen(const Eigen::MatrixXd& measurements, const float& cur_len, const int& index, Eigen::MatrixXd& input_m)
{
	float target_len = (cur_len / measurements.coeff(index, 0))*(input_m.coeff(index, 0));
	return target_len;
}

/*!
*@brief  计算拉普拉斯权值矩阵
*@param[out] 拉普拉斯稀疏矩阵
*@param[in]  const SurfaceMesh & mesh  待求拉普拉斯矩阵的原始网格
*@param[in]  Eigen::SparseMatrix<float> & L  拉普拉斯矩阵，是个大型稀疏矩阵
*@return     void
*/
void FitMeasure::CaculateLaplacianCotMatrix(const SurfaceMesh& mesh, Eigen::SparseMatrix<double> & L, std::vector<Tri>& triplets)
{
	//std::vector<Tri> triplets;
	const int p_num = mesh.n_vertices();
	L.resize(3 * p_num, 3 * p_num);
	for (auto fit : mesh.faces())
	{
		vec3 p[3];
		float cot[3];
		int id[3];
		auto vf = mesh.vertices(fit);
		for (int i = 0; i < 3; ++i, ++vf)
		{
			p[i] = mesh.position(*vf);
			id[i] = (*vf).idx();
		}
		SetTriplets(triplets, cot, p, id);
	}
	L.setFromTriplets(triplets.begin(), triplets.end());
	std::cout << triplets.size() << std::endl;
	std::cout << "Calc Laplacian Done !" << std::endl;
}

void FitMeasure::SetTriplets(std::vector<Tri>& triplets, float cot[3], vec3 p[3], int id[3])
{
	triplets.reserve(7);
	for (int i = 0; i < 3; ++i)
	{
		int j = (i + 1) % 3, k = (j + 1) % 3;
		cot[i] = dot(p[j] - p[i], p[k] - p[i]) / norm(cross(p[j] - p[i], p[k] - p[i]));
		for (int l = 0; l < 3; ++l)
		{
			triplets.emplace_back(Tri(3 * id[j] + l, 3 * id[k] + l, -0.5 * cot[i]));
			triplets.emplace_back(Tri(3 * id[k] + l, 3 * id[j] + l, -0.5 * cot[i]));
			triplets.emplace_back(Tri(3 * id[i] + l, 3 * id[i] + l, 0.5*(cot[(i + 1) % 3] + cot[(i + 2) % 3])));
		}
	}
}

/*!
*@brief  保存每个尺寸边的数量
*@param[out]
*@param[in]  std::vector<std::vector<int>> & point_idx
*@return     void
*/
void FitMeasure::SaveEdge(std::vector<std::vector<int>>& point_idx)
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
void FitMeasure::ConstructCoefficientMatrixBottom(Eigen::SparseMatrix<double>& A, std::vector<std::vector<int>>& point_idx, const Eigen::Matrix3Xd &vertices, const Eigen::MatrixXd& measurements,
	Eigen::VectorXd& b, Eigen::MatrixXd& input_m, std::vector<Tri>& triplets)
{
	//std::vector<Tri> triplets;
	num_verts = vertices.cols();
	num_measure = point_idx.size();	//18
	SaveEdge(point_idx);
	num_edge_all = std::accumulate(edge.begin(), edge.end(), 0);	//所有边的数量1246
	A.resize(3 * num_edge_all, 3 * num_verts);
	SetTriplets(triplets, vertices, b, input_m, point_idx, measurements);
	//std::cout << "b:" << b.rows() << std::endl;
	//A.setFromTriplets(triplets.begin(), triplets.end());
}

void FitMeasure::SetTriplets(std::vector<Tri>& triplets, const Eigen::Matrix3Xd& vertices, Eigen::VectorXd& b, Eigen::MatrixXd& input_m, std::vector<std::vector<int>>& point_idx,
	const Eigen::MatrixXd& measurements)
{
	triplets.reserve(7);
	b.setConstant(3 * num_edge_all, 0);
	int row = 0;
	for (int i = 0; i < num_measure; ++i)
	{
		if (i == 17)
			std::cout << "a" << std::endl;
		for (int j = 0; j < edge[i]; ++j)
		{
			int edge_0;
			int edge_1;
			//i==0身高信息
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
			//std::cout << edge_len << std::endl;
			Eigen::Vector3d auxd = (edge_01 / edge_len) * CalcTargetLen(measurements, edge_len, i, input_m);
			for (int k = 0; k < 3; ++k)
			{
				triplets.push_back(Tri(row + j * 3 + k + 3 * num_verts, edge_0 * 3 + k, -1));
				triplets.push_back(Tri(row + j * 3 + k + 3 * num_verts, edge_1 * 3 + k, 1));
				b(row + j * 3 + k) = auxd[k];
				//std::cout << row + j * 3 + k << " " << edge_0 * 3 + k << " " << auxd[k] << std::endl;
			}
		}
		row += edge[i] * 3;
	}
}

/*!
*@brief  	将顶点矩阵大小3*V转成3V*1
*@param[out]
*@param[in]  Eigen::SparseVector<double> & v
*@param[in]  const Eigen::Matrix3Xd & vertices
*@return     void
*/
void FitMeasure::Mat2Vec(Eigen::SparseVector<double>& v, const Eigen::Matrix3Xd& vertices)
{
	for (int i = 0; i < num_verts; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			v.insert(3 * i + j) = vertices.coeff(0, i);
			v.insert(3 * i + j) = vertices.coeff(1, i);
			v.insert(3 * i + j) = vertices.coeff(2, i);
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
void FitMeasure::ConstructB(Eigen::SparseMatrix<double>& b1, Eigen::VectorXd& b2)
{
	Eigen::VectorXd b(3 * num_verts + 3 * num_edge_all);
	for (int i = 0; i < 3 * num_verts; ++i)
		b(i) = b1.coeff(i, 0);

	for (int j = 0; j < 3 * num_edge_all; ++j)
		b(j + 3 * num_verts) = b2(j);
	ShowMessage(string("FillB"));
}

/*!
*@brief  拟合尺寸
*@param[out]
*@param[in]  Eigen::Matrix3Xd & res_verts
*@param[in]  Eigen::SparseMatrix<double> & C
*@param[in]  Eigen::SparseMatrix<double> & L
*@param[in]  const Eigen::Matrix3Xd & vertices
*@param[in]  Eigen::VectorXd & b2
*@param[in]  std::vector<std::vector<int>> & point_idx
*@return     void
*/
void FitMeasure::FitMeasurements(Eigen::Matrix3Xd& res_verts, Eigen::SparseMatrix<double>& C, Eigen::SparseMatrix<double>& L, const Eigen::Matrix3Xd& vertices, Eigen::VectorXd& b2,
	std::vector<std::vector<int>>& point_idx, std::vector<Tri>& triplets)
{

	Eigen::SparseVector<double> v(3 * num_verts);
	Mat2Vec(v, vertices);

	Eigen::SparseMatrix<double> b1 = L * v;
	ShowMessage(string("b1 = L * v"));

	Eigen::VectorXd b(3 * num_verts + 3 * num_edge_all);
	ConstructB(b1, b2);

	Eigen::SparseMatrix<double> A;
	ConstructCoefficientMatrix(A, L, C,triplets);

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	auto AT = A.transpose();
	ShowMessage(string("AT"));
	std::cout << A.nonZeros() << std::endl;

	Eigen::SparseMatrix<double> ATA = AT * A;
	std::cout << ATA.nonZeros() << std::endl;
	solver.compute(ATA);
	ShowMessage(string("ATA"));

	if (solver.info() != Eigen::Success)
		ShowMessage(string(">Solve Failed"));

	Eigen::VectorXd vecV = solver.solve(A.transpose() * b);
	ShowMessage(string(">Solve Success"));

	res_verts = Eigen::Map<Eigen::Matrix3Xd>(vecV.data(), 3, vertices.cols());
}

void FitMeasure::ShowMessage(const string& msg)
{
	std::cout << msg.c_str() << std::endl;
}

void FitMeasure::ConstructCoefficientMatrix(Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& L, const Eigen::SparseMatrix<double>& C, 
	std::vector<Tri>& triplets)
{
	Eigen::SparseMatrix<double> ttt(3, 3);
	std::vector<Tri> trippp1;
	std::vector<Tri> trippp2;
	trippp1.emplace_back(Tri(1, 1, 2.2));
	ttt.setFromTriplets(trippp1.begin(), trippp1.end());
	ttt.conservativeResize(6, 3);
	trippp2.emplace_back(Tri(5, 2, 33));
	ttt.setFromTriplets(trippp2.begin(), trippp2.end());
	std::cout << "ttt nonezero : " << ttt.nonZeros() << std::endl;

	//std::vector<Tri> triplets;
	//triplets.reserve(7);
	//A = L;
	A.resize(3 * num_verts + 3 * num_edge_all, 3 * num_verts);
	std::cout << "Ainner: " << A.innerSize() << "   " << "Aouter" << A.outerSize() << std::endl;

	//int cnt = 0;
	//for (int i = 3 * num_verts; i < 3 * num_verts + 3 * num_edge_all; ++i)
	//	for (int j = 0; j < 3 * num_verts; ++j)
	//	{
	//		triplets.push_back(Tri(i, j, C.coeff(i - 3 * num_verts, j)));
	//	}
	//for (int i = 0; i < C.outerSize(); ++i)
	//{
	//	for (Eigen::SparseMatrix<double>::InnerIterator it(C, i); it; ++it)
	//	{
	//		triplets.emplace_back(Tri(it.row() + 3 * num_verts, it.col(), it.value()));
	//		//std::cout << it.row() + 3 * num_verts << "  " << it.col() << std::endl;
	//	}
	//}
	A.setFromTriplets(triplets.begin(), triplets.end());



	std::cout << "C nonezero : " << C.nonZeros() << std::endl;
	std::cout << "L nonezero : " << L.nonZeros() << std::endl;
	std::cout << "A nonezero : " << A.nonZeros() << std::endl;
	ShowMessage(string("Construct A"));
}