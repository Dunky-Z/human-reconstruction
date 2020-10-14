#pragma once
#include <Eigen/Dense>

#include "pmp/SurfaceMesh.h"
#include "measure.h"
#include "reshaper.h"
#include "utils.h"

using namespace Eigen;
using namespace pmp;
using namespace std;

typedef Eigen::Triplet<float> Tri;

void CalcEuclideanGradient(Eigen::VectorXd& gradient, Matrix3Xd& vertices);
void CalcGeodesicGradient(Eigen::VectorXd& gradient, Matrix3Xd vertices, Eigen::MatrixXd& measurements);
float CalcTargetLen(Eigen::MatrixXd& measurements, const float& cur_len, const int& index);
void CalcEnergy(double& energy, Eigen::Matrix3Xd& vertices);
void CaculateLaplacianCotMatrix(const SurfaceMesh& mesh, Eigen::SparseMatrix<double> & L);
