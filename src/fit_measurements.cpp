#include<iomanip>

#include "fit_measurements.h"

using namespace Eigen;
using namespace pmp;
using namespace std;


const int M_NUM = 15;
const double step = 1;
Eigen::Matrix3Xd verts;
Eigen::Matrix3Xi faces;
Eigen::MatrixXd input_m(M_NUM, 1);
Eigen::VectorXd gradient;

Eigen::MatrixXd measurements;
std::vector<std::vector<int>> point_idx;
std::vector<std::vector<std::vector<double>>> control_points;


int main()
{
	Reshaper reshaper;
	Measure measure;
	
}

void CalcEnergy(double& enyg, Eigen::Matrix3Xd& vertices)
{
	for (size_t i_ = 1, i = 1; i_ < point_idx.size(); ++i_)
	{
		//剔除三个围长
		if (i_ == 2 || i_ == 3 || i_ == 4)
		{
			continue;
		}
		size_t n = point_idx[i_].size();
		for (size_t j = 0; j < n - 1; ++j)
		{
			int id1 = point_idx[i_][j], id2 = point_idx[i_][j + 1];
			pmp::vec3 p1; p1[0] = vertices.coeff(0, id1); p1[1] = vertices.coeff(1, id1); p1[2] = vertices.coeff(2, id1);
			pmp::vec3 p2; p2[0] = vertices.coeff(0, id2); p2[1] = vertices.coeff(1, id2); p2[2] = vertices.coeff(2, id2);
			double cur_len = distance(p1, p2);
			double p = std::pow(cur_len, 2);
			double l = std::pow(CalcTargetLen(measurements, cur_len, i), 2);
			enyg = enyg + std::pow(p - l, 2);
			std::cout << setprecision(15);
			//std::cout << std::fixed << p << "  " << l << "  " << enyg << std::endl;
		}
		i++;
	}
}

/*!
*@brief  计算欧式距离的梯度
*@param[out]
*@param[in]  Eigen::MatrixXd & gradient
*@param[in]  SurfaceMesh & mesh
*@param[in]  std::vector<int> point_idx
*@param[in]  Eigen::MatrixXd & input_measure
*@return     void
*/
//void CalcEuclideanGradient(Eigen::VectorXd& gradient, Eigen::Matrix3Xd& vertices)
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
void CalcGeodesicGradient(Eigen::VectorXd& gradient, Matrix3Xd vertices, Eigen::MatrixXd& measurements)
{
	pmp::vec3 grad;
	//从1开始因为poin_idx[0]是欧式距离
	for (size_t i_ = 1, i = 1; i_ < point_idx.size(); ++i_)
	{
		if (i_ == 2 || i_ == 3 || i_ == 4)
		{
			continue;
		}
		size_t n = point_idx[i_].size();
		for (size_t j = 0; j < n - 1; ++j)
		{
			int id1 = point_idx[i_][j], id2 = point_idx[i_][j + 1];
			//std::cout << id1 << "  " << id2 << std::endl;
			pmp::vec3 p1; p1[0] = vertices.coeff(0, id1); p1[1] = vertices.coeff(1, id1); p1[2] = vertices.coeff(2, id1);
			pmp::vec3 p2; p2[0] = vertices.coeff(0, id2); p2[1] = vertices.coeff(1, id2); p2[2] = vertices.coeff(2, id2);
			float cur_len = distance(p1, p2);
			grad = 4 * (std::pow(cur_len, 2) - std::pow(CalcTargetLen(measurements, cur_len, i), 2))*(p1 - p2);
			for (size_t k = 0; k < 3; ++k)
			{
				gradient(id1 * 3 + k) += grad[k];
				gradient(id2 * 3 + k) += -grad[k];
			}
		}
		++i;
	}
}

/*!
*@brief  计算目标边长
*@param[out]
*@param[in]  Eigen::MatrixXd measurements
*@param[in]  const float cur_len  当前网格上这条边长
*@param[in]  int index  第index个尺寸
*@return     float
*/
float CalcTargetLen(Eigen::MatrixXd& measurements, const float& cur_len, const int& index)
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
void CaculateLaplacianCotMatrix(const SurfaceMesh& mesh, Eigen::SparseMatrix<double> & L)
{
	std::vector<Tri> tripletlist;
	tripletlist.reserve(20);
	const int p_num = mesh.n_vertices();
	L.resize(p_num, p_num);
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
		for (int i = 0; i < 3; ++i)
		{
			int j = (i + 1) % 3, k = (j + 1) % 3;
			cot[i] = dot(p[j] - p[i], p[k] - p[i]) / norm(cross(p[j] - p[i], p[k] - p[i]));

			tripletlist.push_back(Tri(id[j], id[k], -0.5 * cot[i]));
			tripletlist.push_back(Tri(id[k], id[j], -0.5 * cot[i]));
		}
		for (int i = 0; i < 3; ++i)
		{
			tripletlist.push_back(Tri(id[i], id[i], 0.5*(cot[(i + 1) % 3] + cot[(i + 2) % 3])));
		}
	}
	L.setFromTriplets(tripletlist.begin(), tripletlist.end());
}

void CaculateCoefficientMatrix(const SurfaceMesh& mesh, Eigen::SparseMatrix<double> & A)
{

}


void Reshaper::FitMeasurements(Eigen::Matrix3Xd& res_verts, std::vector<std::vector<int>> point_idx, const Eigen::Matrix3Xd &vertices, const Eigen::MatrixXd measurements)
{
	const int num_measure = point_idx.size() - 1;	//去除第一个体重的信息
	//保存每个尺寸的边数量
	std::vector<int> edge;
	for (auto& num_v : point_idx)
	{
		if (num_v.size() > 2)
		{
			edge.push_back(num_v.size());
		}
		else
		{
			//edge.push_back(1);
			continue;
		}
	}
	//所有边的数量
	int num_edge_all = std::accumulate(edge.begin(), edge.end(), 0);

	typedef Eigen::Triplet<double> Tri;
	std::vector<Tri> triplets;
	triplets.reserve(6 * num_edge_all);
	Eigen::VectorXd b;
	b.setConstant(3 * num_edge_all, 0);
	Eigen::SparseMatrix<double> A;
	A.resize(3 * num_edge_all, 3 * vertices.cols());

	int row = 0;
	for (int i = 0; i < num_measure; ++i)
	{

		for (int j = 0; j < edge[i]; ++j)
		{
			const int edge_0 = point_idx[i + 1][j % edge[i]];
			const int edge_1 = point_idx[i + 1][(j + 1) % edge[i]];
			const Eigen::Vector3d v0 = vertices.col(edge_0);
			const Eigen::Vector3d v1 = vertices.col(edge_1);
			const Eigen::Vector3d edge_01 = v1 - v0;
			const double edge_len = edge_01.norm();
			const Eigen::Vector3d auxd = (edge_01 / edge_len) * (measurements(i + 1, 0) / edge[i]);

			for (int k = 0; k < 3; ++k)
			{
				triplets.push_back(Tri(row + j * 3 + k, edge_0 * 3 + k, -1));
				triplets.push_back(Tri(row + j * 3 + k, edge_1 * 3 + k, 1));
				b(row + j * 3 + k) = auxd[k];
				std::cout << row + j * 3 + k << " " << edge_0 * 3 + k << " " << auxd[k] << std::endl;
			}
		}
		row += edge[i] * 3;
	}

	A.setFromTriplets(triplets.begin(), triplets.end());
	auto AT = A.transpose();
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	//solver.compute(A * AT);
	solver.compute(A.transpose() * A);
	if (solver.info() != Eigen::Success)
		std::cout << "solve failed !" << std::endl;
	Eigen::VectorXd vecV = solver.solve(A.transpose() * b);
	//Eigen::VectorXd vecV = AT *  solver.solve(b);
	res_verts = Eigen::Map<Eigen::Matrix3Xd>(vecV.data(), 3, vertices.cols());
}