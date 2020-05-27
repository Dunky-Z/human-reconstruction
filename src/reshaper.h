#pragma once
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <unsupported/Eigen/KroneckerProduct>

#include "utils.h"
#include "mesh_io.h"
#include "binary_io.h"


class Reshaper
{
public:
	Reshaper();
	~Reshaper();
	void SaveBinControlPoint(std::vector<std::vector<std::vector<double>>> &control_points);
	void SaveBinVerts(const char *filename, const std::string &path, const std::vector<std::string> &files);
	void SaveBinFaces(const char *filename, const std::string &path, const std::vector<std::string> &files);
	void SaveVertFacetInBin();
private:

};
