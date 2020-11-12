#pragma once

#include <vector>
#include <iostream>
#include "pmp/SurfaceMesh.h"
#include <fstream>
#include <Eigen/Cholesky>
#include <windows.h>


#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;
using namespace pmp;
using namespace Eigen;

typedef Eigen::Triplet<double> T;



class LaplaceDeformation
{
public:
	// L: 原始laplace矩阵（方阵），A: 加入了锚点（point_num + anchors, point_num）大小， T：转置
	Eigen::SparseMatrix<double> A, L, AT, ATA, b, ATb, new_vertice;
	std::vector<int>  move_anchor_idx{ 5456 }, fix_anchor_idx{961,118,11410};
	std::vector<Point> move_anchor_coord;							// 形变后移动锚点位置
	Eigen::SparseMatrix<double> AdjacentVertices, VerticesDegree;	//顶点邻接矩阵;顶点度数矩阵：顶点所相邻顶点的数量

	LaplaceDeformation();
	virtual ~LaplaceDeformation() = default;
	void InitializeMesh();
	void AllpyLaplaceDeformation();
	//void BuildAdjacentMatrix(const Surface_mesh &mesh);
	void Compute_CotMatrix(const SurfaceMesh &mesh, std::vector<T>& triplets);
	void BuildATtimesAMatrix(const SurfaceMesh & mesh, std::vector<T>& triplets);
	void BuildATtimesbMatrix(const SurfaceMesh & mesh);
	void SetNewcord(SurfaceMesh & mesh);
	const std::string BIN_DATA_PATH = "../data/";
	const std::string DATASET_PATH = "../data/";
};