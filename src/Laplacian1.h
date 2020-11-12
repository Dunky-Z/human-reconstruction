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
	// L: ԭʼlaplace���󣨷��󣩣�A: ������ê�㣨point_num + anchors, point_num����С�� T��ת��
	Eigen::SparseMatrix<double> A, L, AT, ATA, b, ATb, new_vertice;
	std::vector<int>  move_anchor_idx{ 5456 }, fix_anchor_idx{961,118,11410};
	std::vector<Point> move_anchor_coord;							// �α���ƶ�ê��λ��
	Eigen::SparseMatrix<double> AdjacentVertices, VerticesDegree;	//�����ڽӾ���;����������󣺶��������ڶ��������

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