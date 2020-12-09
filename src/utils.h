#pragma once

#include <io.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <Eigen/Dense>

#define VERTS 12500
#define FACES 25000
#define EPS 1e-6
#define BASIS_NUM 10
#define EPSILON 1e-6

//const std::string DATASET_PATH = "../data/";
const std::string DATASET_PATH = "D:/ITabc/ITabc/objDataSet/dataset/male_tmp/";
const std::string BIN_DATA_PATH = "../data/";
//const std::string BIN_DATA_PATH = "./data/train/male/";

template<class Matrix>
void printShape(Matrix &M) {
	cout << M.rows() << " " << M.cols() << endl;
	return;
}

template<class Matrix, class ...Args>
void printShape(Matrix &M, Args... rest) {
	printShape(M);
	printShape(rest...);
}

std::vector<std::string> GetFiles(const std::string &cate_dir);