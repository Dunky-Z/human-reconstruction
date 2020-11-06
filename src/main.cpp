#include <time.h>
#include <fstream>
#include <iostream>
#include <Eigen/Dense> 


#include "utils.h"
#include "mesh_io.h"
#include "measure.h"
#include "reshaper.h"
#include "fit_measurements.h"

int main()
{
	clock_t t = clock();

	std::vector<std::vector<std::vector<double>>> control_points;
	Reshaper reshaper;
	Measure measure;
	FitMeasure fit;
	SurfaceMesh mesh;

	mesh.read((DATASET_PATH + "1_.obj").c_str());
	reshaper.SaveBinControlPoint(control_points);
	//输入尺寸
	Eigen::MatrixXd input_m(18, 1);
	input_m << 1795.61, 460.47, 1212.81, 1098.78, 1134.35, 890.41, 823.41, 419.05, 824.58, 1126.35, 1199.55, 1336.46, 649.92, 623.889, 204.25, 1313.27, 442.89, 726.47;

	//保存模型顶点，面片信息
	//reshaper.SaveVertFacetInBin();

	Eigen::MatrixXd verts; //shape(3V, num_models)
	Eigen::Matrix3Xi facets;

	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "vertex").c_str(), verts);
	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "facets").c_str(), facets);

	//测量并保存尺寸信息
	Eigen::MatrixXd measurelist;
	//measure.ConvertMeasure(verts, facets, control_points, measurelist);

	Eigen::VectorXd one_measure;//取出一个模型的尺寸信息
	fit.RecoverMeasure(measurelist, one_measure);

	Eigen::Matrix3Xd res_verts;
	Eigen::MatrixXd one_verts = verts.col(0);//取出一个模型的顶点信息
	one_verts.resize(3, verts.rows() / 3);
	//保存边信息
	std::vector<std::vector<int>> point_idx;
	reshaper.SaveBinEdge(control_points, point_idx);

	std::vector<Tri> triplets_A;
	Eigen::SparseMatrix<double> A;
	Eigen::VectorXd b;
	Eigen::SparseMatrix<double>  L;

	fit.CaculateLaplacianCotMatrix(mesh, L, triplets_A);
	fit.ConstructCoefficientMatrixBottom(point_idx, one_verts, one_measure, b, input_m, triplets_A);
	fit.ConstructCoefficientMatrix(A, triplets_A);
	fit.FitMeasurements(res_verts, L, one_verts, b, point_idx, A);
	meshio::SaveObj((BIN_DATA_PATH + "res.obj").c_str(), res_verts, facets);
	cout << res_verts.leftCols(20) << endl;

	cout << "Main spend : " << (double)(clock() - t) / CLOCKS_PER_SEC << "seconds!" << endl;
	getchar();
	return 0;
}