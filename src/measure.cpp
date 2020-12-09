
#include <iostream>

#include "measure.h"

using namespace std;


Measure::Measure() { }

Measure::~Measure() { }


/*!
*@brief  计算一个模型的所有尺寸
*@param[out]
*@param[in]  const std::vector<std::vector<std::vector<double>>> & control_points  控制点信息
*@param[in]  const Eigen::Matrix3Xd & vertices  一个模型的所有顶点信息
*@param[in]  const Eigen::Matrix3Xi & facets  三角面片的信息
*@return     Eigen::MatrixXd
*/Eigen::MatrixXd Measure::CalcMeasure(
	const std::vector<std::vector<std::vector<double>>>& control_points, 
	const Eigen::Matrix3Xd & vertices,
	const Eigen::Matrix3Xi &facets)
{
	Eigen::MatrixXd measure_list;
	measure_list.resize(1, M_NUM);
	double vol = 0.0;
	double kHumanbodyIntensity = 1026.0;
	Eigen::Vector3d v0, v1, v2;
	//cout << F.cols() << endl;
	for (int i = 0; i < facets.cols(); ++i)
	{
		v0 = vertices.col(facets(0, i));
		v1 = vertices.col(facets(1, i));
		v2 = vertices.col(facets(2, i));
		//calc area
		vol += v0.cross(v1).dot(v2);
	}
	vol = abs(vol) / 6.0;


	double weight = kHumanbodyIntensity * vol;
	weight = pow(weight, 1.0 / 3.0) * 1000;
	measure_list(0) = weight;
	int idx = 1;
	//measures[i][j]  第i行数据， j列数据为顶点序号
	for (auto & measures : control_points)
	{
		double length = 0.0;
		Eigen::Vector3d p1, p2;
		p2 = vertices.col(measures[0][1]);

		for (int i = 0; i < measures.size(); ++i)
		{
			p1 = p2;
			if (measures[i][0] == 1)
				p2 = vertices.col(measures[i][1]);
			else if (measures[i][0] == 2)
			{
				p2 = vertices.col(measures[i][1])*measures[i][3] + vertices.col(measures[i][2])*measures[i][4];
			}
			else
			{
				p2 = vertices.col(measures[i][1])*measures[i][4] + vertices.col(measures[i][2])*measures[i][5] + vertices.col(measures[i][3])*measures[i][6];
			}
			length += sqrt((p1 - p2).array().pow(2).sum());
		}
		measure_list(idx++) = length * 1000;
	}
	//shape : 19,1
	return measure_list.transpose();

	//cout << measure_list.cols() << endl;;
}

/*!
*@brief  计算方差
*@param[out]
*@param[in]  const Eigen::MatrixXd & x
*@param[in]  const double average
*@return     double
*/double Measure::CalcVariance(const Eigen::MatrixXd &x, const double average)
{
	double sum = 0.0;
	int len = x.cols();
	for (int i = 0; i < len; ++i)
	{
		sum += pow(x(0, i) - average, 2);
	}
	return sum / len;
}
double Measure::CalcStd(const Eigen::MatrixXd &x, const double average)
{
	double variance = CalcVariance(x, average);
	return sqrt(variance);
}

/*!
*@brief  计算所有模型的尺寸，并以二进制保存
*@param[out]
*@param[in]  const Eigen::MatrixXd & all_vertices  所有模型的顶点信息
*@param[in]  const Eigen::Matrix3Xi & facets  三角面片信息
*@param[in]  const std::vector<std::vector<std::vector<double>>> & control_points  控制点的信息
*@param[in]  Eigen::MatrixXd & measure_lists  测量数据矩阵 shape(19, num_model)
*@return     void
*/void Measure::ConvertMeasure(const Eigen::MatrixXd & all_vertices, const Eigen::Matrix3Xi &facets,
	const std::vector<std::vector<std::vector<double>>>& control_points, Eigen::MatrixXd &measure_lists)
{
	cout << "Start convert measure..." << endl;

	clock_t t = clock();
	Measure measure;
	Eigen::MatrixXd measure_list;
	Eigen::MatrixXd verts;
	measure_lists.resize(19, all_vertices.cols());
	for (int i = 0; i < all_vertices.cols(); ++i)
	{
		cout << i << endl;
		verts = all_vertices.col(i);
		verts.resize(3, 12500);
		measure_list = measure.CalcMeasure(control_points, verts, facets);
		measure_lists.col(i) = measure_list;
	}
	Eigen::MatrixXd mean_measure, std_measure;
	mean_measure.resize(19, 1);
	std_measure.resize(19, 1);

	for (int i = 0; i < 19; ++i)
	{
		mean_measure(i, 0) = measure_lists.row(i).mean();
		std_measure(i, 0) = CalcStd(measure_lists.row(i), mean_measure(i, 0));
	}
	for (int i = 0; i < all_vertices.cols(); ++i)
	{
		measure_lists.col(i) = measure_lists.col(i) - mean_measure;
		measure_lists.col(i) = measure_lists.col(i).cwiseQuotient(std_measure);
	}
	binaryio::WriteMatrixBinaryToFile((BIN_DATA_PATH + "measure_list").c_str(), measure_lists);
	binaryio::WriteMatrixBinaryToFile((BIN_DATA_PATH + "mean_measure").c_str(), mean_measure);
	binaryio::WriteMatrixBinaryToFile((BIN_DATA_PATH + "std_measure").c_str(), std_measure);

	cout << "Calculate measure spend: " << (double)(clock() - t) / CLOCKS_PER_SEC << "seconds." << endl;
}