////#include "Laplacian1.h"
////
////int main(int argc, char **argv)
////{
////	double start = GetTickCount();  //开始时间
////	LaplaceDeformation Deform;
////	Deform.AllpyLaplaceDeformation();
////	double finish = GetTickCount();   //结束时间
////	double t = finish - start;
////	cout << t << endl; //输出时间
////	return 0;
////}
//#include "solver.h"
//
//
//void ShowMessage(const string& msg)
//{
//	std::cout << msg.c_str() << std::endl;
//}
//
//
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
//void ConstructCoefficientMatrix(
//	std::vector<std::vector<int>>& point_idx,
//	const Eigen::Matrix3Xd& vertices,
//	const Eigen::MatrixXd& one_measure,
//	const Eigen::MatrixXd& input_m)
//{
//	//std::vector<Tri> triplets;
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
//void ConstructCoefficientMatrixBottom()
//{
//	C.resize(num_edge_all, num_verts);
//	C.setFromTriplets(triplets_C.begin(), triplets_C.end());
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
//	}
//	for (int i = 0; i < 3; ++i)
//	{
//		triplets_A.push_back({ id[i], id[i], 0.5*(cot[(i + 1) % 3] + cot[(i + 2) % 3]) });
//	}
//}
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
//void ConstructB(
//	const Eigen::SparseMatrix<double>& b1,
//	const Eigen::SparseMatrix<double>& b2)
//{
//	b.resize(num_verts + num_edge_all, 3);
//	for (int i = 0; i < num_verts; ++i)
//	{
//		for (int j = 0; j < 3; ++j)
//		{
//			b.insert(i, j) = b1.coeff(i, j);
//		}
//	}
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
//void solve_one(Eigen::MatrixXd& vertices_one)
//{
//	alglib::real_1d_array x;
//	x.attach_to_ptr(num_verts, vertices_one.data());
//	Eigen::VectorXd scale;
//	scale.resize(num_verts);
//	scale.setOnes();
//	real_1d_array s;
//	s.setcontent(num_verts, scale.data());
//	double epsx = 0.000001;
//	ae_int_t maxits = 2;//最大迭代次数
//	//>设置优化对象，调整一些参数
//	minnlcstate state;
//	minnlccreate(x, state);
//	minnlcsetcond(state, epsx, maxits);
//	minnlcsetscale(state, s);
//	double rho = 1000.0;//惩罚系数
//	ae_int_t outerits = 1;//迭代次数
//	//>选择增广拉格朗日法
//	minnlcsetalgoaul(state, rho, outerits);
//	//>设置约束
//	minnlcsetnlc(state, num_edge_all, 0);
//	minnlcoptguardsmoothness(state);
//	minnlcoptguardgradient(state, 0.001);
//	minnlcreport rep;
//	real_1d_array x1;
//	//>优化求解
//	alglib::minnlcoptimize(state, nlcfunc1_jac);
//	minnlcresults(state, x1, rep);
//	printf("%s\n", x1.tostring(2).c_str());
//	//>检查是否报错
//	optguardreport ogrep;
//	minnlcoptguardresults(state, ogrep);
//	printf("%s\n", ogrep.badgradsuspected ? "true" : "false"); // EXPECTED: false
//	printf("%s\n", ogrep.nonc0suspected ? "true" : "false"); // EXPECTED: false
//	printf("%s\n", ogrep.nonc1suspected ? "true" : "false"); // EXPECTED: false
//}
//
//void solv_three()
//{
//	for (int i = 0; i < 3; ++i)
//	{
//		Eigen::MatrixXd vertices_one = vertices.row(i);
//		b_col = i;
//		solve_one(vertices_one);
//	}
//}
//
//void  nlcfunc1_jac(const real_1d_array& x, real_1d_array& fi, real_2d_array& jac, void* ptr)
//{
//	ShowMessage("Calc Jacobian");
//	alglib::real_1d_array temp(x);
//	Eigen::MatrixXd V0 = Eigen::Map<Eigen::MatrixXd>(temp.getcontent(), num_verts, 1);
//	Eigen::MatrixXd F_up = L * V0;
//	Eigen::MatrixXd F_down = C * V0;
//
//	for (int i = 0; i < num_verts; ++i)
//	{
//		auto t = F_up.coeff(i, 0) - b_up.coeff(i, b_col);
//		fi[0] += t;
//	}
//	for (int i = 1; i <= num_edge_all; ++i)
//	{
//		auto t = F_down.coeff(i - 1, 0) - b_down.coeff(i - 1, b_col);
//		fi[i] = t;
//	}
//	auto LT = L.transpose();
//	auto CT = C.transpose();
//	Eigen::SparseMatrix<double> jacob_up = 2 * LT * (L * V0 - b_up.col(b_col)+ b_up.col(b_col));
//	for (int i = 0; i < num_verts; ++i)
//	{
//		jac[0][i] = jacob_up.coeff(i, 0);
//	}
//	Eigen::SparseMatrix<double> jacob_down = 2 * CT * (C * V0 - b_down.col(b_col)+ b_down.col(b_col));
//	for (int i = 1; i <= num_edge_all; ++i)
//	{
//		for (int j = 0; j < num_verts; ++j)
//		{
//			jac[i][j] = jacob_down.coeff(j, 0);
//		}
//	}
//}
//
//int main(int argc, char **argv)
//{
//	FitMeasure	fit;
//	SurfaceMesh mesh;
//	Measure		measure;
//	Reshaper	reshaper;
//	//输入尺寸
//	Eigen::MatrixXd input_m(18, 1);
//	input_m << 1695.61, 460.47, 1312.81, 1098.78, 1134.35, 890.41, 823.41, 419.05,
//		824.58, 1126.35, 1299.55, 1336.46, 649.92, 623.889, 204.25, 1313.27, 442.89, 726.47;
//
//	reshaper.SaveBinControlPoint(control_points);
//	mesh.read((DATASET_PATH + "1_.obj").c_str());
//	meshio::ReadObj((DATASET_PATH + "1_.obj").c_str(), vertices, faces);
//	Eigen::MatrixXd one_measure;
//	one_measure.resize(19, 1);
//	one_measure = measure.CalcMeasure(control_points, vertices, faces);
//	std::vector<std::vector<int>> point_idx;
//	reshaper.SaveBinEdge(control_points, point_idx);
//	reshaper.SaveBinControlPoint(control_points);
//	CaculateLaplacianCotMatrix(mesh);
//	ConstructCoefficientMatrix(point_idx, vertices, one_measure, input_m);
//	ConstructCoefficientMatrixBottom();
//	Eigen::SparseMatrix<double> v(num_verts, 3);
//	Mat2Vec(v, vertices);
//
//	Eigen::SparseMatrix<double> b1 = L * v;
//	b_up = b1;
//	ShowMessage(string("b1 = L * v"));
//	ConstructB(b1, b_down);
//
//	Eigen::MatrixXd vertices_one = vertices.row(0);
//	//solv_three();
//	solve_one(vertices_one);
//	return 0;
//}
