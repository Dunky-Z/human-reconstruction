#pragma once
#pragma once
#include <time.h>
#include <string>
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseQR>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>	

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../alglib/cpp/src/optimization.h"
#include "utils.h"
#include "measure.h"
#include "reshaper.h"

#include "pmp/SurfaceMesh.h"

using namespace pmp;
using namespace std;
using namespace Eigen;
using namespace alglib;

typedef Eigen::Triplet<double> Tri;



void CalcEuclideanGradient(Eigen::VectorXd& gradient, Matrix3Xd& vertices);
void CalcGeodesicGradient(Eigen::VectorXd& gradient, Matrix3Xd& vertices, Eigen::MatrixXd& measurements, Eigen::MatrixXd& input_m);
void CalcEnergy(double& energy, Eigen::Matrix3Xd& vertices);
void CaculateLaplacianCotMatrix_Test(
	SurfaceMesh& mesh,
	Eigen::SparseMatrix<double> & L,
	std::vector<Tri>& triplets,
	Eigen::SparseMatrix<double>& b);
void ConstructCoefficientMatrixBottom_Test(
	SurfaceMesh& mesh,
	std::vector<std::vector<int>>& point_idx,
	const Eigen::MatrixXd& measurements,
	Eigen::SparseMatrix<double>& b,
	const Eigen::MatrixXd& input_m,
	std::vector<Tri>& triplets);
void FitMeasurements_Test(
	SurfaceMesh& mesh,
	Eigen::Matrix3Xd& res_verts,
	const Eigen::SparseMatrix<double>& L,
	const Eigen::Matrix3Xd& vertices,
	Eigen::SparseMatrix<double>&	b2,
	std::vector<std::vector<int>>& point_idx,
	Eigen::SparseMatrix<double> A);
void Mat2Vec_Test(
	Eigen::SparseMatrix<double>& v,
	const std::vector<Point>& Vertice);
float CalcTargetLen(
	const Eigen::MatrixXd& measurements,
	const float& cur_len,
	const int& index,
	const Eigen::MatrixXd& input_m);
void CaculateLaplacianCotMatrix(
	SurfaceMesh& mesh,
	std::vector<Tri>& triplets_A);
void ConstructCoefficientMatrix(
	std::vector<std::vector<int>>& point_idx,
	Eigen::Matrix3Xd& vertices,
	const Eigen::MatrixXd& measurements,
	const Eigen::MatrixXd& input_m,
	std::vector<Tri>& triplets_A,
	std::vector<Tri>& triplets_C);
void FitMeasurements(
	SurfaceMesh& mesh,
	Eigen::Matrix3Xd& vertices);
void SaveEdge(std::vector<std::vector<int>>& point_idx);
void SetTriplets(
	vec3 p[3],
	int id[3],
	std::vector<Tri>& triplets_A);
void SetTriplets(
	Eigen::Matrix3Xd& vertices,
	const Eigen::MatrixXd& input_m,
	const std::vector<std::vector<int>>& point_idx,
	const Eigen::MatrixXd& measurements,
	std::vector<Tri>& triplets_A,
	std::vector<Tri>& triplets_C);
void Mat2Vec(Eigen::SparseMatrix<double>& v,
	Eigen::Matrix3Xd& vertices);
void ConstructB(
	Eigen::SparseMatrix<double>& b1,
	Eigen::SparseMatrix<double>& b2);
void ShowMessage(
	const string& msg);
void ConstructCoefficientMatrixBottom(
	std::vector<Tri>& triplets_C);
void RecoverMeasure(
	Eigen::MatrixXd& measurelist,
	Eigen::MatrixXd& one_measure);
void SaveObj(
	SurfaceMesh& mesh,
	Eigen::SparseMatrix<double>& new_vertice);
void SaveObj(
	SurfaceMesh& mesh,
	Eigen::MatrixX3d& new_vertice);
void SaveObj(
	SurfaceMesh& mesh,
	Eigen::Matrix3Xd& new_vertice);
Eigen::SparseMatrix<double> CalcGradient(
	Eigen::MatrixXd& V0,
	Eigen::SparseMatrix<double>& LL,
	Eigen::SparseMatrix<double>& bb);
Eigen::MatrixXd SolveOneDimension(
	Eigen::MatrixXd& vertices_one);
Eigen::MatrixX3d AULSolver(
	Eigen::Matrix3Xd& vertices);
void SetGrad(
	Eigen::MatrixXd& grad_t,
	real_1d_array &grad);
Eigen::MatrixXd LBFGS(
	SurfaceMesh& mesh,
	Eigen::Matrix3Xd& vertices,
	std::vector<Tri>& triplets_A,
	std::vector<Tri>& triplets_C);
void array2mat(
	real_1d_array& x,
	Eigen::MatrixXd& res);
void array2mat(
	real_1d_array& x,
	Eigen::Matrix3Xd& res);
Eigen::Matrix3Xd LeastSquare(
	SurfaceMesh& mesh,
	Eigen::Matrix3Xd& vertices,
	std::vector<Tri>& triplets_A,
	std::vector<Tri>& triplets_C);
Eigen::MatrixXd Solver(
	SurfaceMesh& mesh,
	Eigen::Matrix3Xd& vertices);

void Mat2Array(
	Eigen::Matrix3Xd& vertices,
	Eigen::MatrixXd& vertices_t);
double func_x(
	Eigen::SparseMatrix<double>& LL,
	Eigen::MatrixXd& Vp,
	Eigen::SparseMatrix<double>& bb);
Eigen::MatrixXd Mat2Array(
	Eigen::MatrixXd& v);
Eigen::MatrixXd Mat2Array(
	Eigen::SparseMatrix<double>& v);
void function1_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr);
