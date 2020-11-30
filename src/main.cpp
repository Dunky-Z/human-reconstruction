//#include <time.h>
//#include <fstream>
//#include <iostream>
//#include <Eigen/Dense> 
//
//
//#include "utils.h"
//#include "mesh_io.h"
//#include "measure.h"
//#include "reshaper.h"
//#include "fit_measurements.h"
//
//int main()
//{
//	clock_t t = clock();
//
//	std::vector<std::vector<std::vector<double>>> control_points;
//
//	FitMeasure	fit;
//	SurfaceMesh mesh;
//	Measure		measure;
//	Reshaper	reshaper;
//
//	fit.FitMeasurements();
//	mesh.read((DATASET_PATH + "1_.obj").c_str());
//	reshaper.SaveBinControlPoint(control_points);
//	//输入尺寸
//	Eigen::MatrixXd input_m(18, 1);
//	input_m << 1695.61, 460.47, 1312.81, 1098.78, 1134.35, 890.41, 823.41, 419.05, 
//		824.58, 1126.35, 1299.55, 1336.46, 649.92, 623.889, 204.25, 1313.27, 442.89, 726.47;
//
//	Eigen::MatrixXd verts; //shape(3V, num_models)
//	Eigen::Matrix3Xi facets;
//
//	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "vertex").c_str(), verts);
//	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "facets").c_str(), facets);
//
//	Eigen::MatrixXd one_verts = verts.col(0);//取出一个模型的顶点信息
//	one_verts.resize(3, verts.rows() / 3);
//
//	Eigen::MatrixXd one_measure;
//	one_measure.resize(19, 1);
//	one_measure = measure.CalcMeasure(control_points, one_verts, facets);
//
//	std::vector<std::vector<int>> point_idx;
//	reshaper.SaveBinEdge(control_points, point_idx);
//
//	fit.FitMeasurements(mesh, one_verts, one_measure, point_idx, input_m);
//
//	cout << "Main spend : " << (double)(clock() - t) / CLOCKS_PER_SEC << "seconds!" << endl;
//	getchar();
//	return 0;
//}