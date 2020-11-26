//#include "Laplacian1.h"
//
//int main(int argc, char **argv)
//{
//	double start = GetTickCount();  //开始时间
//	LaplaceDeformation Deform;
//	Deform.AllpyLaplaceDeformation();
//	double finish = GetTickCount();   //结束时间
//	double t = finish - start;
//	cout << t << endl; //输出时间
//	return 0;
//}

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fit_measurements.h"
#include "../alglib/cpp/src/optimization.h"

using namespace pmp;
using namespace std;
using namespace Eigen;
using namespace alglib;

typedef Eigen::Triplet<double> Tri;

Eigen::Matrix3Xd vertices;
Eigen::Matrix3Xi faces;
Eigen::SparseMatrix<double> A;
Eigen::SparseMatrix<double> L;
Eigen::SparseMatrix<double> b;
Eigen::SparseMatrix<double> b_down;
int num_verts;
int num_measure;
int num_edge_all;
std::vector<int> edge;

std::vector<Tri> triplets_A;
std::vector<std::vector<std::vector<double>>> control_points;


void ShowMessage(const string& msg)
{
	std::cout << msg.c_str() << std::endl;
}

void ConstructCoefficientMatrix()
{
	A.resize(num_verts + num_edge_all, num_verts);
	A.setFromTriplets(triplets_A.begin(), triplets_A.end());
	ShowMessage(string("Construct A"));
}

float CalcTargetLen(
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

void SaveEdge(
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

void ConstructCoefficientMatrixBottom(
	std::vector<std::vector<int>>& point_idx,
	const Eigen::Matrix3Xd &vertices,
	const Eigen::MatrixXd& one_measure,
	const Eigen::MatrixXd& input_m)
{
	//std::vector<Tri> triplets;
	num_verts = vertices.cols();
	num_measure = point_idx.size();	//18
	SaveEdge(point_idx);
	num_edge_all = std::accumulate(edge.begin(), edge.end(), 0);	//所有边的数量1246
	b_down.resize(num_edge_all, 3);
	SetTriplets(vertices, input_m, point_idx, one_measure);
}

void SetTriplets(
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
			for (int k = 0; k < 3; ++k)
			{
				b_down.insert(row + j, k) = auxd[k];
			}
			//std::cout << "row: " << row + j * 3 + k << " " << "col: " << edge_0 * 3 + k << " " << auxd[k] << std::endl;
		}
		row += edge[i];
	}
	std::cout << row << std::endl;
}

void SetTriplets(
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
void CaculateLaplacianCotMatrix(
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

void ConstructB(
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

void Mat2Vec(
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

int main(int argc, char **argv)
{
	FitMeasure	fit;
	SurfaceMesh mesh;
	Measure		measure;
	Reshaper	reshaper;


	//输入尺寸
	Eigen::MatrixXd input_m(18, 1);
	input_m << 1695.61, 460.47, 1312.81, 1098.78, 1134.35, 890.41, 823.41, 419.05,
		824.58, 1126.35, 1299.55, 1336.46, 649.92, 623.889, 204.25, 1313.27, 442.89, 726.47;

	reshaper.SaveBinControlPoint(control_points);
	mesh.read((DATASET_PATH + "1_.obj").c_str());
	meshio::ReadObj((DATASET_PATH + "1_.obj").c_str(), vertices, faces);
	Eigen::MatrixXd one_measure;
	one_measure.resize(19, 1);
	one_measure = measure.CalcMeasure(control_points, vertices, faces);

	std::vector<std::vector<int>> point_idx;
	reshaper.SaveBinEdge(control_points, point_idx);

	reshaper.SaveBinControlPoint(control_points);
	CaculateLaplacianCotMatrix(mesh);
	ConstructCoefficientMatrixBottom(point_idx, vertices, one_measure, input_m);
	Eigen::SparseMatrix<double> v(num_verts, 3);
	Mat2Vec(v, vertices);

	Eigen::SparseMatrix<double> b1 = L * v;
	ShowMessage(string("b1 = L * v"));
	ConstructB(b1, b_down);


	alglib::real_1d_array x;
	x.attach_to_ptr(num_verts * 3, vertices.data());
	//
	// This example demonstrates minimization of
	//
	//     f(x0,x1) = -x0+x1
	//
	// subject to nonlinear equality constraint
	//
	//    x0^2 + x1^2 - 1 = 0
	//
	real_1d_array s = "[1,1]";
	double epsx = 0.000001;
	ae_int_t maxits = 0;
	minnlcstate state;

	//
	// Create optimizer object and tune its settings:
	// * epsx=0.000001  stopping condition for inner iterations
	// * s=[1,1]        all variables have unit scale
	//
	minnlccreate(x, state);
	minnlcsetcond(state, epsx, maxits);
	minnlcsetscale(state, s);

	//
	// Choose one of the nonlinear programming solvers supported by minnlc
	// optimizer:
	// * SLP - successive linear programming NLP solver
	// * AUL - augmented Lagrangian NLP solver
	//
	// Different solvers have different properties:
	// * SLP is the most robust solver provided by ALGLIB: it can solve both
	//   convex and nonconvex optimization problems, it respects box and
	//   linear constraints (after you find feasible point it won't move away
	//   from the feasible area) and tries to respect nonlinear constraints
	//   as much as possible. It also usually needs less function evaluations
	//   to converge than AUL.
	//   However, it solves LP subproblems at each iterations which adds
	//   significant overhead to its running time. Sometimes it can be as much
	//   as 7x times slower than AUL.
	// * AUL solver is less robust than SLP - it can violate box and linear
	//   constraints at any moment, and it is intended for convex optimization
	//   problems (although in many cases it can deal with nonconvex ones too).
	//   Also, unlike SLP it needs some tuning (penalty factor and number of
	//   outer iterations).
	//   However, it is often much faster than the current version of SLP.
	//
	// In the code below we set solver to be AUL but then override it with SLP,
	// so the effective choice is to use SLP. We recommend you to use SLP at
	// least for early prototyping stages.
	//
	// You can comment out line with SLP if you want to solve your problem with
	// AUL solver.
	//
	double rho = 1000.0;
	ae_int_t outerits = 5;
	minnlcsetalgoaul(state, rho, outerits);
	//minnlcsetalgoslp(state);

	//
	// Set constraints:
	//
	// Nonlinear constraints are tricky - you can not "pack" general
	// nonlinear function into double precision array. That's why
	// minnlcsetnlc() does not accept constraints itself - only constraint
	// counts are passed: first parameter is number of equality constraints,
	// second one is number of inequality constraints.
	//
	// As for constraining functions - these functions are passed as part
	// of problem Jacobian (see below).
	//
	// NOTE: MinNLC optimizer supports arbitrary combination of boundary, general
	//       linear and general nonlinear constraints. This example does not
	//       show how to work with general linear constraints, but you can
	//       easily find it in documentation on minnlcsetbc() and
	//       minnlcsetlc() functions.
	//
	minnlcsetnlc(state, 1, 0);

	//
	// Activate OptGuard integrity checking.
	//
	// OptGuard monitor helps to catch common coding and problem statement
	// issues, like:
	// * discontinuity of the target/constraints (C0 continuity violation)
	// * nonsmoothness of the target/constraints (C1 continuity violation)
	// * erroneous analytic Jacobian, i.e. one inconsistent with actual
	//   change in the target/constraints
	//
	// OptGuard is essential for early prototyping stages because such
	// problems often result in premature termination of the optimizer
	// which is really hard to distinguish from the correct termination.
	//
	// IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
	//            DIFFERENTIATION, THUS DO NOT USE IT IN PRODUCTION CODE!
	//
	//            Other OptGuard checks add moderate overhead, but anyway
	//            it is better to turn them off when they are not needed.
	//
	minnlcoptguardsmoothness(state);
	minnlcoptguardgradient(state, 0.001);

	//
	// Optimize and test results.
	//
	// Optimizer object accepts vector function and its Jacobian, with first
	// component (Jacobian row) being target function, and next components
	// (Jacobian rows) being nonlinear equality and inequality constraints.
	//
	// So, our vector function has form
	//
	//     {f0,f1} = { -x0+x1 , x0^2+x1^2-1 }
	//
	// with Jacobian
	//
	//         [  -1    +1  ]
	//     J = [            ]
	//         [ 2*x0  2*x1 ]
	//
	// with f0 being target function, f1 being constraining function. Number
	// of equality/inequality constraints is specified by minnlcsetnlc(),
	// with equality ones always being first, inequality ones being last.
	//
	minnlcreport rep;
	real_1d_array x1;
	alglib::minnlcoptimize(state, nlcfunc1_jac);
	minnlcresults(state, x1, rep);
	printf("%s\n", x1.tostring(2).c_str()); // EXPECTED: [0.70710,-0.70710]

	//
	// Check that OptGuard did not report errors
	//
	// NOTE: want to test OptGuard? Try breaking the Jacobian - say, add
	//       1.0 to some of its components.
	//
	optguardreport ogrep;
	minnlcoptguardresults(state, ogrep);
	printf("%s\n", ogrep.badgradsuspected ? "true" : "false"); // EXPECTED: false
	printf("%s\n", ogrep.nonc0suspected ? "true" : "false"); // EXPECTED: false
	printf("%s\n", ogrep.nonc1suspected ? "true" : "false"); // EXPECTED: false
	return 0;
}


void  nlcfunc1_jac(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr)
{
	//
	// this callback calculates
	//
	//     f0(x0,x1) = -x0+x1
	//     f1(x0,x1) = x0^2+x1^2-1
	//
	// and Jacobian matrix J = [dfi/dxj]
	//
	alglib::real_1d_array temp(x);
	Eigen::Map<Matrix3Xd> V0(temp.getcontent(), 3, num_verts);
	auto F = A * V0;

	for (int i = 0; i < 3 * (num_verts + num_edge_all); ++i)
	{
		fi[i] = F.coeff(i / 3, i % 3) - b.coeff(i / 3, i % 3);
	}
	auto AT = A.transpose();
	auto ATA = AT * A;
	Eigen::SparseMatrix<double> jacob = 2*AT*(A*V0 - b);
	for (int i = 0; i < (num_verts+num_edge_all); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			jac[i][j] = jacob.coeff(i, j);
		}
	}
}