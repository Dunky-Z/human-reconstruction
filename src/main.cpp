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
	//����ߴ�
	Eigen::MatrixXd input_m(18, 1);
	input_m << 1695.61, 460.47, 1312.81, 1098.78, 1134.35, 890.41, 823.41, 419.05, 824.58, 1126.35, 1299.55, 1336.46, 649.92, 623.889, 204.25, 1313.27, 442.89, 726.47;

	//����ģ�Ͷ��㣬��Ƭ��Ϣ
	//reshaper.SaveVertFacetInBin();

	Eigen::MatrixXd verts; //shape(3V, num_models)
	Eigen::Matrix3Xi facets;

	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "vertex").c_str(), verts);
	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "facets").c_str(), facets);

	//����������ߴ���Ϣ
	Eigen::MatrixXd measurelist;
	//measure.ConvertMeasure(verts, facets, control_points, measurelist);

	Eigen::VectorXd one_measure;//ȡ��һ��ģ�͵ĳߴ���Ϣ
	fit.RecoverMeasure(measurelist, one_measure);

	Eigen::Matrix3Xd res_verts;
	Eigen::MatrixXd one_verts = verts.col(0);//ȡ��һ��ģ�͵Ķ�����Ϣ
	one_verts.resize(3, verts.rows() / 3);
	//�������Ϣ
	std::vector<std::vector<int>> point_idx;
	reshaper.SaveBinEdge(control_points, point_idx);

	//fit.CaculateLaplacianCotMatrix_Test(mesh, L, triplets_A, b);
	fit.CaculateLaplacianCotMatrix(mesh);
	fit.ConstructCoefficientMatrixBottom(point_idx, one_verts, one_measure, input_m);
	//fit.ConstructCoefficientMatrixBottom_Test(mesh,point_idx, one_verts, b, input_m, triplets_A);
	fit.ConstructCoefficientMatrix();
	fit.FitMeasurements(mesh,res_verts, one_verts, point_idx);
	meshio::SaveObj((BIN_DATA_PATH + "res.obj").c_str(), res_verts, facets);
	cout << res_verts.leftCols(20) << endl;

	cout << "Main spend : " << (double)(clock() - t) / CLOCKS_PER_SEC << "seconds!" << endl;
	getchar();
	return 0;
}