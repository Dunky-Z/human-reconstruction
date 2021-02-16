#include <time.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <Eigen/Dense> 
#include "vtk.h"
#include "vtk_reader.h"
#include "reshaper.h"
#include "measure.h"
#include "pmp/SurfaceMesh.h"
#include "utils.h"


using namespace pmp;
using namespace std;

std::vector<std::vector<int>> point_idx;
std::vector<std::vector<std::vector<double>>> control_points;



void LengthError()
{
	std::vector<std::vector<std::vector<double>>> control_points;
}

template<typename OS>
void writeLines(OS& os, const Eigen::Matrix3Xd &V, std::vector<int>& p_idx, bool circle)
{
	long seg = V.cols() - 1;
	std::vector<int> line_idx;
	if (circle) seg++;
	line_idx.reserve(seg);
	for (long i = 0; i < seg; i++)
	{
		line_idx.emplace_back(i);
		line_idx.emplace_back((i + 1) % V.cols());
	}
	line2vtk(os, V.data(), V.cols(), line_idx.data(), seg);
}

Eigen::Matrix3Xd SaveCurveVerts(const Eigen::Matrix3Xd& verts, std::vector<int>& p_idx)
{
	Eigen::Matrix3Xd curve_verts;
	curve_verts.resize(3, p_idx.size());
	int i = 0;
	for (auto& idx : p_idx)
	{
		curve_verts.coeffRef(0, i) = verts.coeff(0, idx);
		curve_verts.coeffRef(1, i) = verts.coeff(1, idx);
		curve_verts.coeffRef(2, i) = verts.coeff(2, idx);
		i++;
	}
	return curve_verts;
}

int main()
{
	SurfaceMesh mesh;
	Measure		measure;
	Reshaper	reshaper;

	mesh.read((BIN_DATA_PATH + "AVE.obj").c_str());

	Eigen::Matrix3Xd verts;
	Eigen::Matrix3Xi facets;

	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "facets").c_str(), facets);
	meshio::ReadObj((BIN_DATA_PATH + "AVE.obj").c_str(), verts, facets);

	reshaper.SaveBinControlPoint(control_points);
	reshaper.SaveBinEdge(control_points, point_idx);

	//std::vector<int> p_idx = point_idx[18];
	for (int i = 0; i < 18; ++i)
	{
		for (auto& aa : point_idx[i])
		{
			cout << aa << endl;
		}
		cout << "-------" << endl;
	}

	Eigen::Matrix3Xd cur_verts = SaveCurveVerts(verts, p_idx);
	cout << cur_verts << endl;
	int seg = p_idx.size() + 1;
	std::ofstream os((BIN_DATA_PATH + "Maximum thigh.vtk").c_str());
	writeLines(os, cur_verts, p_idx, 1);
	return 0;
}

