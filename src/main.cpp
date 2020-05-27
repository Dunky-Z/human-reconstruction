#include <time.h>
#include <fstream>
#include <iostream>
#include <Eigen/Dense> 

#include "utils.h"
#include "measure.h"
#include "reshaper.h"

int main()
{
	clock_t t = clock();

	std::vector<std::vector<std::vector<double>>> control_points;
	Reshaper reshaper;
	Measure measure;
	reshaper.SaveBinControlPoint(control_points);
	//保存模型顶点，面片信息
	//reshaper.SaveVertFacetInBin();

	Eigen::MatrixXd verts;
	Eigen::Matrix3Xi facets;

	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "vertex").c_str(), verts);
	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "facets").c_str(), facets);

	//测量并保存尺寸信息
	Eigen::MatrixXd measurelist;
	measure.ConvertMeasure(verts, facets, control_points, measurelist);
	cout << "Main spend : " << (double)(clock() - t) / CLOCKS_PER_SEC << "seconds!" << endl;
	getchar();
	return 0;
}