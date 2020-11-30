#pragma once

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
Eigen::SparseMatrix<double> C;
Eigen::SparseMatrix<double> b;
Eigen::SparseMatrix<double> b_down;
Eigen::SparseMatrix<double> b_up;
int b_col;
int num_verts;
int num_measure;
int num_edge_all;

std::vector<int> edge;

std::vector<Tri> triplets_A;
std::vector<Tri> triplets_C;
std::vector<std::vector<std::vector<double>>> control_points;

void ShowMessage(const string& msg);
void ConstructCoefficientMatrix(
	std::vector<std::vector<int>>& point_idx,
	const Eigen::Matrix3Xd& vertices,
	const Eigen::MatrixXd& one_measure,
	const Eigen::MatrixXd& input_m);
float CalcTargetLen(
	const Eigen::MatrixXd& one_measure,
	const float& cur_len,
	const int& index,
	const Eigen::MatrixXd& input_m);
void SaveEdge(
	std::vector<std::vector<int>>& point_idx);
void ConstructCoefficientMatrixBottom();
void SetTriplets(
	const Eigen::Matrix3Xd& vertices,
	const Eigen::MatrixXd& input_m,
	const std::vector<std::vector<int>>& point_idx,
	const Eigen::MatrixXd& one_measure);
void SetTriplets(
	vec3 p[3],
	int id[3]);
void CaculateLaplacianCotMatrix(
	const SurfaceMesh& mesh);
void ConstructB(
	const Eigen::SparseMatrix<double>& b1,
	const Eigen::SparseMatrix<double>& b2);
void Mat2Vec(
	Eigen::SparseMatrix<double>& v,
	const Eigen::Matrix3Xd& vertices);
void  nlcfunc1_jac(
	const real_1d_array &x,
	real_1d_array &fi,
	real_2d_array &jac, 
	void *ptr);