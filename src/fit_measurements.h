#pragma once
#include <time.h>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include<Eigen/SparseQR>
#include<Eigen/SparseLU>
#include<Eigen/SparseCholesky>	

#include "pmp/SurfaceMesh.h"
#include "measure.h"
#include "reshaper.h"
#include "utils.h"

using namespace Eigen;
using namespace pmp;
using namespace std;

typedef Eigen::Triplet<float> Tri;

void CalcEuclideanGradient(Eigen::VectorXd& gradient, Matrix3Xd& vertices);
void CalcGeodesicGradient(Eigen::VectorXd& gradient, Matrix3Xd vertices, Eigen::MatrixXd& measurements, Eigen::MatrixXd& input_m);
float CalcTargetLen(const Eigen::MatrixXd& measurements, const float& cur_len, const int& index, Eigen::MatrixXd& input_m);
void CalcEnergy(double& energy, Eigen::Matrix3Xd& vertices);
void CaculateLaplacianCotMatrix(const SurfaceMesh& mesh, Eigen::SparseMatrix<double> & L);
void CaculateCoefficientMatrix(Eigen::SparseMatrix<double>& A, std::vector<std::vector<int>>& point_idx, const Eigen::Matrix3Xd &vertices, const Eigen::MatrixXd& measurements,
	Eigen::VectorXd& b, Eigen::MatrixXd& input_m);
void FitMeasurements(Eigen::Matrix3Xd& res_verts, Eigen::SparseMatrix<double>& C, Eigen::SparseMatrix<double>& L, const Eigen::Matrix3Xd &vertices, Eigen::VectorXd& b2,
	std::vector<std::vector<int>>& point_idx);