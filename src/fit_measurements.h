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

class FitMeasure
{
public:
	FitMeasure();
	~FitMeasure();

	void CalcEuclideanGradient(Eigen::VectorXd& gradient, Matrix3Xd& vertices);
	void CalcGeodesicGradient(Eigen::VectorXd& gradient, Matrix3Xd& vertices, Eigen::MatrixXd& measurements, Eigen::MatrixXd& input_m);
	void CalcEnergy(double& energy, Eigen::Matrix3Xd& vertices);
	void FitMeasure::CaculateLaplacianCotMatrix_Test(
		SurfaceMesh& mesh,
		Eigen::SparseMatrix<double> & L,
		std::vector<Tri>& triplets,
		Eigen::SparseMatrix<double>& b);
	void FitMeasure::ConstructCoefficientMatrixBottom_Test(
		SurfaceMesh& mesh,
		std::vector<std::vector<int>>& point_idx,
		const Eigen::MatrixXd& measurements,
		Eigen::SparseMatrix<double>& b,
		const Eigen::MatrixXd& input_m,
		std::vector<Tri>& triplets);
	void FitMeasure::FitMeasurements_Test(
		SurfaceMesh& mesh,
		Eigen::Matrix3Xd& res_verts,
		const Eigen::SparseMatrix<double>& L,
		const Eigen::Matrix3Xd& vertices,
		Eigen::SparseMatrix<double>&	b2,
		std::vector<std::vector<int>>& point_idx,
		Eigen::SparseMatrix<double> A);
	void FitMeasure::Mat2Vec_Test(
		Eigen::SparseMatrix<double>& v,
		const std::vector<Point>& Vertice);
	float CalcTargetLen(
		const Eigen::MatrixXd& measurements, 
		const float& cur_len, 
		const int& index, 
		const Eigen::MatrixXd& input_m);
	void CaculateLaplacianCotMatrix(
		const SurfaceMesh& mesh);
	void ConstructCoefficientMatrix(
		std::vector<std::vector<int>>& point_idx,
		const Eigen::Matrix3Xd& vertices,
		const Eigen::MatrixXd& measurements,
		const Eigen::MatrixXd& input_m);
	void FitMeasurements(
		SurfaceMesh& mesh,
		const Eigen::Matrix3Xd& vertices,
		const Eigen::MatrixXd& measurements,
		std::vector<std::vector<int>>& point_idx,
		const Eigen::MatrixXd& input_m);
	void SaveEdge(std::vector<std::vector<int>>& point_idx);
	void SetTriplets(
		vec3 p[3],
		int id[3]);
	void SetTriplets(
		const Eigen::Matrix3Xd& vertices,
		const Eigen::MatrixXd& input_m,
		const std::vector<std::vector<int>>& point_idx,
		const Eigen::MatrixXd& measurements);
	void Mat2Vec(Eigen::SparseMatrix<double>& v,
		const Eigen::Matrix3Xd& vertices);
	void ConstructB(
		Eigen::SparseMatrix<double>& b1,
		Eigen::SparseMatrix<double>& b2);
	void ShowMessage(const string& msg);
	void ConstructCoefficientMatrixBottom();
	void FitMeasure::RecoverMeasure(
		Eigen::MatrixXd& measurelist,
		Eigen::MatrixXd& one_measure);
	void FitMeasure::SaveObj(
		SurfaceMesh& mesh,
		Eigen::SparseMatrix<double>& new_vertice);
	void FitMeasure::SaveObj(
		SurfaceMesh& mesh,
		Eigen::MatrixX3d& new_vertice);
	Eigen::SparseMatrix<double> FitMeasure::CalcGradient(
		Eigen::MatrixXd& V0);
	Eigen::MatrixXd FitMeasure::SolveOneDimension(
		Eigen::MatrixXd& vertices_one);
	Eigen::MatrixX3d FitMeasure::AULSolver(
		const Eigen::Matrix3Xd& vertices);
	void SetGrad(Eigen::MatrixXd& grad_t, real_1d_array &grad);
	Eigen::MatrixXd FitMeasure::LBFGS(
		Eigen::MatrixXd& vertices_one);
	void FitMeasure::array2mat(
		real_1d_array& x, 
		Eigen::MatrixXd& res);

	friend void function1_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr);

protected:
private:
	int M_NUM;
	const double step = 1;
	const double eps = 1e-6;
	double delta = 0.01;
	int b_col = 0;
	int num_measure;
	int num_edge_all;
	int num_verts;
	double beta = 1;
	Eigen::MatrixXd alpha;
	std::vector<Tri> triplets_A;
	std::vector<Tri> triplets_C;
	std::vector<int> edge;//size = 18
	std::vector<Point> Vertice;
	Eigen::Matrix3Xd verts;
	Eigen::Matrix3Xi faces;
	Eigen::VectorXd gradient;
	Eigen::SparseMatrix<double> A;
	Eigen::SparseMatrix<double> L;
	Eigen::SparseMatrix<double> C;
	Eigen::SparseMatrix<double> b;
	Eigen::SparseMatrix<double> b_up;
	Eigen::SparseMatrix<double> b_down;
	//std::vector<std::vector<int>> point_idx;//[m][n]m个尺寸，每个尺寸n个点，每个点对应的顶点下标
	std::vector<std::vector<std::vector<double>>> control_points;

};