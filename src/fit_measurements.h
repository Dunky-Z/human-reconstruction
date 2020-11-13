#pragma once
#include <time.h>
#include <string>
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

typedef Eigen::Triplet<double> Tri;

class FitMeasure
{
public:
	FitMeasure();
	~FitMeasure();

	void CalcEuclideanGradient(Eigen::VectorXd& gradient, Matrix3Xd& vertices);
	void CalcGeodesicGradient(Eigen::VectorXd& gradient, Matrix3Xd& vertices, Eigen::MatrixXd& measurements, Eigen::MatrixXd& input_m);
	float CalcTargetLen(const Eigen::MatrixXd& measurements, const float& cur_len, const int& index, const Eigen::MatrixXd& input_m);
	void CalcEnergy(double& energy, Eigen::Matrix3Xd& vertices);
	void CaculateLaplacianCotMatrix(
		const SurfaceMesh& mesh);
	void ConstructCoefficientMatrixBottom(
		std::vector<std::vector<int>>& point_idx, 
		const Eigen::Matrix3Xd &vertices,
		const Eigen::MatrixXd& measurements,
		const Eigen::MatrixXd& input_m);
	void FitMeasurements(
		SurfaceMesh& mesh,
		Eigen::Matrix3Xd& res_verts,
		const Eigen::Matrix3Xd& vertices, 
		std::vector<std::vector<int>>& point_idx);
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
	void ConstructCoefficientMatrix();
	void FitMeasure::RecoverMeasure(Eigen::MatrixXd& measurelist, Eigen::VectorXd& one_measure);
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
		Eigen::SparseMatrix<double> & b2,
		std::vector<std::vector<int>>& point_idx,
		Eigen::SparseMatrix<double> A);
	void FitMeasure::Mat2Vec_Test(
		Eigen::SparseMatrix<double>& v,
		const std::vector<Point>& Vertice);
protected:
private:
	int M_NUM;
	const double step = 1;
	int num_measure;
	int num_edge_all;
	int num_verts;
	std::vector<int> edge;//size = 18
	Eigen::Matrix3Xd verts;
	Eigen::Matrix3Xi faces;
	Eigen::VectorXd gradient;
	//std::vector<std::vector<int>> point_idx;//[m][n]m���ߴ磬ÿ���ߴ�n���㣬ÿ�����Ӧ�Ķ����±�
	std::vector<std::vector<std::vector<double>>> control_points;
	std::vector<Point> Vertice;
	Eigen::SparseMatrix<double> A;
	std::vector<Tri> triplets_A;
	Eigen::SparseMatrix<double> L;
	Eigen::SparseMatrix<double>& b_down;
	Eigen::SparseMatrix<double> b;

};