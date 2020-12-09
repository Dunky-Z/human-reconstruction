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

	FitMeasure	fit;
	SurfaceMesh mesh;
	Measure		measure;
	Reshaper	reshaper;

	mesh.read((DATASET_PATH + "R1.obj").c_str());
	reshaper.SaveBinControlPoint(control_points);
	//输入尺寸
	Eigen::MatrixXd input_m(18, 1);
	//input_m << 1990.61, 440.47, 1190.76, 1143.59, 1134.35, 942.00, 957.41, 413.05,
	//	868.58, 1106.35, 1226.55, 1386.46, 708.92, 404.889, 201.25, 1424.27, 447.89, 800.47;

	//input_m << 1927.511930, 417.407466, 1080.334960, 1020.748780, 1080.516651, 935.276469, 900.714513, 426.964179, 874.499698, 985.313128, 1181.424510,
	//	1297.044574, 686.606327, 401.353873, 189.919517, 1402.991225, 426.907171, 706.888987;

	//input_m << 1780, 400.2769, 1071.95308, 1140, 1130, 940, 910, 404.81295, 860, 1100, 1225.90881, 1306.96954, 700,
	//	400, 200, 1323.72555, 440, 800;

	input_m << 1795.61, 460.47, 1212.81, 1098.78, 1134.35, 890.41, 823.41, 419.05, 824.58, 1126.354,
		1199.55, 1336.46, 649.92, 423.889, 204.25, 1313.27, 442.89, 726.47;

	Eigen::Matrix3Xd verts; //shape(3V, num_models)
	Eigen::Matrix3Xi facets;

	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "facets").c_str(), facets);
	meshio::ReadObj((DATASET_PATH + "R3.obj").c_str(), verts, facets);

	Eigen::MatrixXd one_verts = verts;//取出一个模型的顶点信息
	one_verts.resize(3, verts.cols());

	Eigen::MatrixXd one_measure;
	one_measure.resize(19, 1);
	one_measure = measure.CalcMeasure(control_points, one_verts, facets);
	std::cout << one_measure << std::endl;
	std::vector<std::vector<int>> point_idx;
	reshaper.SaveBinEdge(control_points, point_idx);

	fit.FitMeasurements(mesh, one_verts, one_measure, point_idx, input_m);
	fit.ErrorAnalysis(measure, control_points, input_m);
	cout << "Main spend : " << (double)(clock() - t) / CLOCKS_PER_SEC << "seconds!" << endl;
	getchar();
	return 0;
}