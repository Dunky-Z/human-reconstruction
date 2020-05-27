#pragma once

#include <algorithm>
#include <Eigen/Dense>

#include "utils.h"
#include "binary_io.h"

using namespace std;
class Measure
{
private:
	static const int M_NUM = 19;

public:

	Measure();
	~Measure();

	double CalcVariance(const Eigen::MatrixXd &x, const double average);
	double CalcStd(const Eigen::MatrixXd &x, const double average);
	Eigen::MatrixXd Measure::CalcMeasure(const std::vector<std::vector<std::vector<double>>>& control_points, const Eigen::Matrix3Xd & vertices,
		const Eigen::Matrix3Xi &facets);
	void Measure::ConvertMeasure(const Eigen::MatrixXd & all_vertices, const Eigen::Matrix3Xi &facets,
		const std::vector<std::vector<std::vector<double>>>& control_points, Eigen::MatrixXd &measure_lists);
};

