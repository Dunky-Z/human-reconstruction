//#include<iomanip>
//
//#include "fit_measure.h"
//
//using namespace pmp;
//using namespace std;
//
//
//
//int M_NUM;
//const double step = 1;
//const double eps = 1e-6;
//double delta = 0.01;
//int b_col = 0;
//int num_measure;
//int num_edge_all;
//int num_verts;
//double beta = 1;
//Eigen::MatrixXd alpha;
//std::vector<Tri> triplets_A;
//std::vector<Tri> triplets_C;
//std::vector<int> edge;//size = 18
//std::vector<Point> Vertice;
//Eigen::Matrix3Xd verts;
//Eigen::Matrix3Xi faces;
//Eigen::VectorXd gradient;
//Eigen::SparseMatrix<double> A;
//Eigen::SparseMatrix<double> L;
//Eigen::SparseMatrix<double> C;
//Eigen::SparseMatrix<double> b;
//Eigen::SparseMatrix<double> b_up;
//Eigen::SparseMatrix<double> b_down;
////std::vector<std::vector<int>> point_idx;//[m][n]m个尺寸，每个尺寸n个点，每个点对应的顶点下标
//std::vector<std::vector<std::vector<double>>> control_points;
//
//
//#include <time.h>
//#include <fstream>
//#include <iostream>
//#include <Eigen/Dense> 
//
//
//#include "utils.h"
//#include "mesh_io.h"
//#include "measure.h"
//#include "reshaper.h"
//
//int main()
//{
//	clock_t t = clock();
//
//	std::vector<std::vector<std::vector<double>>> control_points;
//
//	SurfaceMesh mesh;
//	Measure		measure;
//	Reshaper	reshaper;
//
//	mesh.read((DATASET_PATH + "1_.obj").c_str());
//	reshaper.SaveBinControlPoint(control_points);
//	//输入尺寸
//	Eigen::MatrixXd input_m(18, 1);
//	input_m << 1695.61, 460.47, 1312.81, 1098.78, 1134.35, 890.41, 823.41, 419.05,
//		824.58, 1126.35, 1299.55, 1336.46, 649.92, 623.889, 204.25, 1313.27, 442.89, 726.47;
//
//	Eigen::MatrixXd verts; //shape(3V, num_models)
//	Eigen::Matrix3Xi facets;
//
//	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "vertex").c_str(), verts);
//	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "facets").c_str(), facets);
//
//	Eigen::MatrixXd one_verts = verts.col(0);//取出一个模型的顶点信息
//	one_verts.resize(3, verts.rows() / 3);
//
//	Eigen::MatrixXd one_measure;
//	one_measure.resize(19, 1);
//	one_measure = measure.CalcMeasure(control_points, one_verts, facets);
//
//	std::vector<std::vector<int>> point_idx;
//	reshaper.SaveBinEdge(control_points, point_idx);
//
//	FitMeasurements(mesh, one_verts, one_measure, point_idx, input_m);
//
//	cout << "Main spend : " << (double)(clock() - t) / CLOCKS_PER_SEC << "seconds!" << endl;
//	getchar();
//	return 0;
//}
//
///*!
//*@brief  计算目标边长
//*@param[out]
//*@param[in]  Eigen::MatrixXd measurements
//*@param[in]  const float cur_len  当前网格上这条边长
//*@param[in]  int index  第index个尺寸
//*@return     float
//*/
//float CalcTargetLen(
//	const Eigen::MatrixXd& one_measure,
//	const float& cur_len,
//	const int& index,
//	const Eigen::MatrixXd& input_m)
//{
//	float m = one_measure.coeff(index + 1, 0) / 1000.0;
//	float cur_len_p = cur_len / m;
//	float target_len_p = input_m.coeff(index, 0) / 1000.0;
//	float target_len = cur_len_p * target_len_p;
//	return target_len;
//}
//
///*!
//*@brief  计算拉普拉斯权值矩阵
//*@param[out] 拉普拉斯稀疏矩阵
//*@param[in]  const SurfaceMesh & mesh  待求拉普拉斯矩阵的原始网格
//*@param[in]  Eigen::SparseMatrix<float> & L  拉普拉斯矩阵，是个大型稀疏矩阵
//*@return     void
//*/
//void CaculateLaplacianCotMatrix(
//	const SurfaceMesh& mesh)
//{
//	//std::vector<Tri> triplets;
//	const int p_num = mesh.n_vertices();
//	L.resize(p_num, p_num);
//	auto fit = mesh.faces_begin();
//	do {
//		auto vf = mesh.vertices(*fit);
//		Point p[3];
//		int id[3];
//		for (int i = 0; i < 3; ++i, ++vf)
//		{
//			p[i] = mesh.position(*vf);
//			id[i] = (*vf).idx();
//		}
//		SetTriplets(p, id);
//	} while (++fit != mesh.faces_end());
//	L.setFromTriplets(triplets_A.begin(), triplets_A.end());
//	std::cout << "Calc Laplacian Done !" << std::endl;
//}
//
//void SetTriplets(
//	vec3 p[3],
//	int id[3])
//{
//	triplets_A.reserve(7);
//	float cot[3];
//	for (int i = 0; i < 3; ++i)
//	{
//		int j = (i + 1) % 3, k = (j + 1) % 3;
//		cot[i] = dot(p[j] - p[i], p[k] - p[i]) / norm(cross(p[j] - p[i], p[k] - p[i]));
//		triplets_A.push_back({ id[j],id[k], -0.5 * cot[i] });
//		triplets_A.push_back({ id[k],id[j], -0.5 * cot[i] });
//
//	}
//	for (int i = 0; i < 3; ++i)
//	{
//		triplets_A.push_back({ id[i], id[i], 0.5*(cot[(i + 1) % 3] + cot[(i + 2) % 3]) });
//	}
//
//}
//
///*!
//*@brief  保存每个尺寸边的数量
//*@param[out]
//*@param[in]  std::vector<std::vector<int>> & point_idx
//*@return     void
//*/
//void SaveEdge(
//	std::vector<std::vector<int>>& point_idx)
//{
//	//保存每个尺寸的边数量
//	//std::vector<int> edge;//size = 18
//	for (auto& num_v : point_idx)
//	{
//		if (num_v.size() > 2)
//		{
//			edge.push_back(num_v.size());
//		}
//		else
//		{
//			edge.push_back(1);
//			//continue;
//		}
//	}
//}
//
//void ConstructCoefficientMatrixBottom()
//{
//	C.resize(num_edge_all, num_verts);
//	C.setFromTriplets(triplets_C.begin(), triplets_C.end());
//}
///*!
//*@brief  构建系数矩阵，A矩阵下半部分
//*@param[out]
//*@param[in]  Eigen::SparseMatrix<double> & A
//*@param[in]  std::vector<std::vector<int>> & point_idx
//*@param[in]  const Eigen::Matrix3Xd & vertices
//*@param[in]  const Eigen::MatrixXd & measurements
//*@param[in]  Eigen::VectorXd & b
//*@param[in]  Eigen::MatrixXd & input_m
//*@return     void
//*/
//
//void ConstructCoefficientMatrix(
//	std::vector<std::vector<int>>& point_idx,
//	const Eigen::Matrix3Xd &vertices,
//	const Eigen::MatrixXd& one_measure,
//	const Eigen::MatrixXd& input_m)
//{
//	num_verts = vertices.cols();
//	num_measure = point_idx.size();	//18
//	SaveEdge(point_idx);
//	num_edge_all = std::accumulate(edge.begin(), edge.end(), 0);	//所有边的数量1246
//	b_down.resize(num_edge_all, 3);
//	SetTriplets(vertices, input_m, point_idx, one_measure);
//	A.resize(num_verts + num_edge_all, num_verts);
//	A.setFromTriplets(triplets_A.begin(), triplets_A.end());
//	ShowMessage(string("Construct A"));
//}
//
//void SetTriplets(
//	const Eigen::Matrix3Xd& vertices,
//	const Eigen::MatrixXd& input_m,
//	const std::vector<std::vector<int>>& point_idx,
//	const Eigen::MatrixXd& one_measure)
//{
//	triplets_A.reserve(7);
//	int row = 0;
//	for (int i = 0; i < num_measure; ++i)
//	{
//		//遍历每个尺寸的每条边
//		for (int j = 0; j < edge[i]; ++j)
//		{
//			int edge_0;
//			int edge_1;
//			//i == 0身高信息
//			if (i == 0)
//			{
//				edge_0 = point_idx[i][j];
//				edge_1 = point_idx[i][j + 1];
//			}
//			else
//			{
//				edge_0 = point_idx[i][j % edge[i]];
//				edge_1 = point_idx[i][(j + 1) % edge[i]];
//			}
//			const Eigen::Vector3d v0 = vertices.col(edge_0);
//			const Eigen::Vector3d v1 = vertices.col(edge_1);
//			const Eigen::Vector3d edge_01 = v1 - v0;
//			const double edge_len = edge_01.norm();
//			double target_len = CalcTargetLen(one_measure, edge_len, i, input_m);
//			//std::cout << edge_len << std::endl;
//			Eigen::Vector3d auxd = (edge_01 / edge_len) * target_len;
//
//			int _row = num_verts + row + j;
//			int _col0 = edge_0;
//			int _col1 = edge_1;
//			triplets_A.emplace_back(Tri(_row, _col0, -1));
//			triplets_A.emplace_back(Tri(_row, _col1, 1));
//			triplets_C.emplace_back(Tri(row + j, _col0, -1));
//			triplets_C.emplace_back(Tri(row + j, _col1, 1));
//			for (int k = 0; k < 3; ++k)
//			{
//				b_down.insert(row + j, k) = auxd[k];
//			}
//			//std::cout << "row: " << row + j * 3 + k << " " << "col: " << edge_0 * 3 + k << " " << auxd[k] << std::endl;
//		}
//		row += edge[i];
//	}
//	std::cout << row << std::endl;
//}
//
///*!
//*@brief  	将顶点矩阵大小3*V转成3V*1
//*@param[out]
//*@param[in]  Eigen::SparseVector<double> & v
//*@param[in]  const Eigen::Matrix3Xd & vertices
//*@return     void
//*/
//void Mat2Vec(
//	Eigen::SparseMatrix<double>& v,
//	const Eigen::Matrix3Xd& vertices)
//{
//	for (int i = 0; i < num_verts; ++i)
//	{
//		for (int j = 0; j < 3; ++j)
//		{
//			v.coeffRef(i, j) = vertices.coeff(j, i);
//		}
//	}
//}
//
///*!
//*@brief  构建矩阵b，即方程右边
//*@param[out]
//*@param[in]  Eigen::SparseMatrix<double> & b1
//*@param[in]  Eigen::VectorXd & b2
//*@return     void
//*/
//void ConstructB(
//	Eigen::SparseMatrix<double>& b1,
//	Eigen::SparseMatrix<double>& b2)
//{
//	b.resize(num_verts + num_edge_all, 3);
//	for (int i = 0; i < num_verts; ++i)
//		for (int j = 0; j < 3; ++j)
//		{
//			b.insert(i, j) = b1.coeff(i, j);
//		}
//
//	for (int i = 0; i < num_edge_all; ++i)
//	{
//		for (int j = 0; j < 3; ++j)
//		{
//			b.insert(num_verts + i, j) = b2.coeff(i, j);
//		}
//	}
//	ShowMessage(string("FillB"));
//}
//
//
//void SaveObj(
//	SurfaceMesh& mesh,
//	Eigen::SparseMatrix<double>& new_vertice)
//{
//	for (auto v : mesh.vertices())
//	{
//		Point& pos = mesh.position(v);
//		pos[0] = new_vertice.coeff(v.idx(), 0);
//		pos[1] = new_vertice.coeff(v.idx(), 1);
//		pos[2] = new_vertice.coeff(v.idx(), 2);
//	}
//	mesh.write((BIN_DATA_PATH + "res.obj").c_str());
//}
//void SaveObj(
//	SurfaceMesh& mesh,
//	Eigen::MatrixX3d& new_vertice)
//{
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
//void ShowMessage(const string& msg)
//{
//	std::cout << msg.c_str() << std::endl;
//}
//
//
//
////将归一化的尺寸恢复为原来大小
//void RecoverMeasure(Eigen::MatrixXd& measurelist, Eigen::MatrixXd& one_measure)
//{
//	Eigen::MatrixXd mean_measure, std_measure;
//	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "mean_measure").c_str(), mean_measure);
//	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "std_measure").c_str(), std_measure);
//	one_measure = one_measure.cwiseProduct(std_measure);
//	one_measure += mean_measure;
//}
//
//Eigen::MatrixX3d AULSolver(
//	const Eigen::Matrix3Xd& vertices)
//{
//	Eigen::MatrixX3d res_vert; res_vert.resize(num_verts, 3);
//	for (int i = 0; i < 3; ++i)
//	{
//		Eigen::MatrixXd vertices_one = vertices.row(i);
//		b_col = i;
//		res_vert.col(i) = SolveOneDimension(vertices_one);
//	}
//	return res_vert;
//}
//
//Eigen::MatrixXd SolveOneDimension(
//	Eigen::MatrixXd& vertices_one)
//{
//	vertices_one.resize(num_verts, 1);
//	auto CT = C.transpose();
//	Eigen::MatrixXd vertices_one_next = vertices_one;
//	Eigen::MatrixXd alpha_next(alpha);
//	alpha.resize(num_edge_all, 1); alpha.setConstant(0.1);
//
//	do
//	{
//		alpha = alpha_next;
//		vertices_one_next = LBFGS(vertices_one);
//		alpha_next = alpha + beta * (C * vertices_one_next - b_down.col(b_col));
//	} while (vertices_one.norm() - vertices_one_next.norm() < eps);
//	return vertices_one_next;
//}
//
//Eigen::SparseMatrix<double> CalcGradient(
//	Eigen::MatrixXd& V0)
//{
//	Eigen::MatrixXd F_up = L * V0;
//	Eigen::MatrixXd F_down = C * V0;
//	auto LT = L.transpose();
//	auto CT = C.transpose();
//	Eigen::SparseMatrix<double> grad_fx = 2 * LT * (L * V0 - b_up.col(b_col) + b_up.col(b_col));
//	return grad_fx;
//}
//
//void SetGrad(Eigen::MatrixXd& grad_t, real_1d_array& grad)
//{
//	for (int i = 0; i < grad_t.size(); ++i)
//	{
//		grad[i] = grad_t.coeff(i, 0);
//	}
//}
//
//Eigen::MatrixXd LBFGS(Eigen::MatrixXd& vertices_one)
//{
//	alglib::real_1d_array x;
//	x.attach_to_ptr(num_verts, vertices_one.data());
//	Eigen::VectorXd scale;
//	scale.resize(num_verts);
//	scale.setOnes();
//	real_1d_array s;
//	s.setcontent(num_verts, scale.data());
//
//	double epsg = 0;
//	double epsf = 0;
//	double epsx = 1e-8;
//	ae_int_t maxits = 0;
//	minlbfgsstate state;
//	minlbfgscreate(1, x, state);
//	minlbfgssetcond(state, epsg, epsf, epsx, maxits);
//	minlbfgssetscale(state, s);
//
//
//	minlbfgsoptguardsmoothness(state);
//	minlbfgsoptguardgradient(state, 0.001);
//	minlbfgsreport rep;
//
//	alglib::minlbfgsoptimize(state, function1_grad);
//	minlbfgsresults(state, x, rep);
//
//	Eigen::MatrixXd res;
//	array2mat(x, res);
//
//	return res;
//
//	//
//	// Check that OptGuard did not report errors
//	//
//	// NOTE: want to test OptGuard? Try breaking the gradient - say, add
//	//       1.0 to some of its components.
//	//
//	optguardreport ogrep;
//	minlbfgsoptguardresults(state, ogrep);
//	printf("%s\n", ogrep.badgradsuspected ? "true" : "false"); // EXPECTED: false
//	printf("%s\n", ogrep.nonc0suspected ? "true" : "false"); // EXPECTED: false
//	printf("%s\n", ogrep.nonc1suspected ? "true" : "false"); // EXPECTED: false
//}
//
//void array2mat(real_1d_array& x, Eigen::MatrixXd& res)
//{
//	res.resize(num_verts, 1);
//	for (int i = 0; i < num_verts; ++i)
//	{
//		res.coeffRef(i, 0) = x[i];
//	}
//}
//
//void function1_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr)
//{
//
//	alglib::real_1d_array temp(x);
//	Eigen::MatrixXd V0 = Eigen::Map<Eigen::MatrixXd>(temp.getcontent(), num_verts, 1);
//	Eigen::MatrixXd F_up = L * V0;
//	Eigen::MatrixXd F_down = C * V0;
//	auto CT = C.transpose();
//	alpha.resize(num_edge_all, 1); alpha.setConstant(0.1);
//
//	for (int i = 0; i < num_verts; ++i)
//	{
//		double t = F_up.coeff(i, 0) - b_up.coeff(i, b_col);
//		func += t;
//	}
//	for (int i = 1; i <= num_edge_all; ++i)
//	{
//		auto t = F_down.coeff(i - 1, 0) - b_down.coeff(i - 1, b_col);
//		func += t;
//	}
//
//	Eigen::MatrixXd grad_t = CalcGradient(V0) + CT * alpha +
//		beta * CT*(C * V0 - b_down.col(b_col));
//	SetGrad(grad_t, grad);
//}
//
//
///*!
//*@brief  拟合尺寸
//*@param[out]
//*@param[in]  Eigen::Matrix3Xd & res_verts
//*@param[in]  Eigen::SparseMatrix<double> & C
//*@param[in]  const Eigen::Matrix3Xd & vertices
//*@param[in]  Eigen::VectorXd & b2
//*@param[in]  std::vector<std::vector<int>> & point_idx
//*@return     void
//*/
//void FitMeasurements(
//	SurfaceMesh& mesh,
//	const Eigen::Matrix3Xd& vertices,
//	const Eigen::MatrixXd& one_measure,
//	std::vector<std::vector<int>>& point_idx,
//	const Eigen::MatrixXd& input_m)
//{
//	CaculateLaplacianCotMatrix(mesh);
//	ConstructCoefficientMatrix(point_idx, vertices, one_measure, input_m);
//	ConstructCoefficientMatrixBottom();
//
//	Eigen::SparseMatrix<double> v(num_verts, 3);
//	Mat2Vec(v, vertices);
//
//	Eigen::SparseMatrix<double> b1 = L * v;
//	b_up = b1;
//	ShowMessage(string("b1 = L * v"));
//	ConstructB(b1, b_down);
//	Eigen::MatrixX3d new_vertice = AULSolver(vertices);
//	SaveObj(mesh, new_vertice);
//	/*Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
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
//	SaveObj(mesh, new_vertice);*/
//}
