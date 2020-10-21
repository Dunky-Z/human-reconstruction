#include <time.h>
#include <fstream>
#include <iostream>
#include <Eigen/Dense> 


#include "utils.h"
#include "measure.h"
#include "reshaper.h"
#include "fit_measurements.h"

int main()
{
	clock_t t = clock();

	std::vector<std::vector<std::vector<double>>> control_points;
	Reshaper reshaper;
	Measure measure;
	SurfaceMesh mesh;
	mesh.read((DATASET_PATH + "1_.obj").c_str());
	reshaper.SaveBinControlPoint(control_points);
	//保存模型顶点，面片信息
	//reshaper.SaveVertFacetInBin();

	Eigen::MatrixXd verts; //shape(3V, num_models)
	Eigen::Matrix3Xi facets;

	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "vertex").c_str(), verts);
	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "facets").c_str(), facets);

	//测量并保存尺寸信息
	Eigen::MatrixXd measurelist;
	//measure.ConvertMeasure(verts, facets, control_points, measurelist);
	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "measure_list").c_str(), measurelist);
	Eigen::VectorXd one_measure = measurelist.col(0);//取出一个模型的尺寸信息
	Eigen::Matrix3Xd res;
	Eigen::MatrixXd one_verts = verts.col(0);//取出一个模型的顶点信息
	one_verts.resize(3, verts.rows() / 3);
	//保存边信息
	std::vector<std::vector<int>> point_idx;
	reshaper.SaveBinEdge(control_points, point_idx);

	Eigen::SparseMatrix<double> C;
	Eigen::VectorXd b;
	CaculateCoefficientMatrix(C, point_idx, one_verts, one_measure, b);
	FitMeasurements(res, point_idx, one_verts, one_measure);

	cout << res.leftCols(20) << endl;


	cout << "Main spend : " << (double)(clock() - t) / CLOCKS_PER_SEC << "seconds!" << endl;
	getchar();
	return 0;
}