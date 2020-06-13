#include "reshaper.h"


Reshaper::Reshaper() {}

Reshaper::~Reshaper() {}


/*!
*@brief  ��ȡ���Ƶ��ļ����������Ƶ㱣�������
*@param[out]
*@param[in]  std::vector<std::vector<std::vector<double>>> & control_points  ������Ƶ������
*@return     void
*/void Reshaper::SaveBinControlPoint(std::vector<std::vector<std::vector<double>>>& control_points)
{
	std::cout << "Begin load cp..." << std::endl;
	std::ifstream is("./data/body_control_points.txt");
	assert(is);

	std::vector<std::vector<double>> templist;
	std::vector<double> tempNode;
	std::string line;
	double t_num;
	while (!is.eof())
	{
		std::getline(is, line);
		if (line.empty())
			continue;

		if (line[0] == '#')
		{
			//cout << "ÿ���ߴ��м����㣺 " << templist.size() << endl;
			if (templist.size() != 0)
			{
				control_points.push_back(templist);
				templist.clear();
			}
		}
		else if (line.find(" ") == std::string::npos)
		{
			continue;
		}
		else
		{
			std::istringstream instream(line);
			//���������Կո�Ϊ�ָ� ����һ������
			while (instream >> t_num)
			{
				//cout << t_num << "  ";
				tempNode.push_back(t_num);
			}
			//cout << "ÿ�����м���Ԫ�� �� "<< tempNode.size() << endl;
			templist.push_back(tempNode);
			tempNode.clear();
		}
	}
	//cout << "ÿ���ߴ��м����㣺 " << templist.size() << endl;
	control_points.push_back(templist);
	std::cout << "Load control_points done!" << std::endl;
}

/*!
*@brief
*@param[out]
*@param[in]  std::vector<std::vector<std::vector<double>>> & control_points  [�ߴ������ÿ���ߴ�����Ŀ��Ƶ㣬ÿ�����Ƶ�]
*@param[in]  const char * filename
*@param[in]  const std::string & path
*@param[in]  const std::vector<std::string> & files
*@return     void
*/
void Reshaper::SaveBinEdge(std::vector<std::vector<std::vector<double>>>& control_points, std::vector<std::vector<int>> &point_idx)
{
	Eigen::Matrix2Xi edge;
	//��Ϊcontrolpoints���������ظ��Ŀ��Ƶ���Ϣ������Ҫȥ�أ�ʣ�����Ķ���Ϳɴ����ߣ������������
	for (auto &measure : control_points)
	{
		std::vector<int> point_idx_;
		for (int i = 0; i < measure.size(); ++i)
		{
			if (measure[i][0] == 1)
			{
				//�ڶ��������в����Ƿ��Ѿ����ڸö��㣬������ھ�����
				std::vector<int>::iterator it = std::find(point_idx_.begin(), point_idx_.end(), measure[i][1]);
				if (it == point_idx_.end())
				{
					point_idx_.push_back(measure[i][1]);
				}
			}
			else if (measure[i][0] == 2)
			{
				auto it1 = std::find(point_idx_.begin(), point_idx_.end(), measure[i][1]);
				if (it1 == point_idx_.end())
				{
					point_idx_.push_back(measure[i][1]);
				}
				auto it2 = std::find(point_idx_.begin(), point_idx_.end(), measure[i][2]);
				if (it2 == point_idx_.end())
				{
					point_idx_.push_back(measure[i][2]);
				}
			}
			else
			{
				std::vector<int>::iterator it1 = std::find(point_idx_.begin(), point_idx_.end(), measure[i][1]);
				if (it1 == point_idx_.end())
				{
					point_idx_.push_back(measure[i][1]);
				}
				std::vector<int>::iterator it2 = std::find(point_idx_.begin(), point_idx_.end(), measure[i][2]);
				if (it2 == point_idx_.end())
				{
					point_idx_.push_back(measure[i][2]);
				}
				std::vector<int>::iterator it3 = std::find(point_idx_.begin(), point_idx_.end(), measure[i][3]);
				if (it3 == point_idx_.end())
				{
					point_idx_.push_back(measure[i][3]);
				}
			}
		}
		point_idx.push_back(point_idx_);
	}
}

void Reshaper::FitOneMeasurements(Eigen::Matrix3Xd & res, std::vector<int> point_idx, const Eigen::Matrix3Xd & vertices, const double measurement)
{
	const int num_v = point_idx.size();
	int num_edge;
	if (num_v > 2)
		num_edge = num_v;
	else
		num_edge = 1;

	Eigen::SparseMatrix<double> A;
	Eigen::VectorXd b;
	typedef Eigen::Triplet<double> Tri;
	std::vector<Tri> triplets;

	A.resize(num_edge, 3 * num_v);
	b.setConstant(num_edge, 0);

	for (int i = 0; i < 3 * num_edge; i = i + 3)
	{
		const int edge_0 = point_idx[i % (3 * num_edge)];
		const int edge_1 = point_idx[(i + 1) % (3 * num_edge)];
		const Eigen::Vector3d& v0 = vertices.col(edge_0);
		const Eigen::Vector3d& v1 = vertices.col(edge_1);
		const Eigen::Vector3d& edge_01 = v1 - v0;
		const double edge_len = edge_01.norm();

		for (int j = 0; j < 3; ++j)
		{
			triplets.push_back(Tri(i + j, edge_0 * 3, -1));
			triplets.push_back(Tri(i + j, edge_1 * 3, 1));
		}
		Eigen::Vector3d one;
		one << 1, 1, 1;
		one.transpose();
		b(i) = one.dot((edge_01 / edge_len)*(measurement / num_edge));
	}




	A.setFromTriplets(triplets.begin(), triplets.end());
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(A.transpose() * A);
	Eigen::VectorXd vecV = solver.solve(A.transpose() * b);
	res = Eigen::Map<Eigen::Matrix3Xd>(vecV.data(), 3, num_v);
}
//
//void Reshaper::FitMeasurements(std::vector<std::vector<int>> point_idx, const Eigen::Matrix3Xd &vertices, const Eigen::MatrixXd measurements)
//{
//	int num_measure = point_idx.size();
//	Eigen::MatrixXd res_verts;
//	for (int i = 0; i < num_measure; ++i)
//	{
//		FitOneMeasurements(point_idx[i], vertices, measurements.col(i));
//	}
//}
/*!
*@brief  �Զ����Ʊ��涥����Ϣ
*@param[out]
*@param[in]  const char * filename  ������ļ���
*@param[in]  const std::string & path  ģ�����ݼ�·��
*@param[in]  const std::vector<std::string> & files  ģ�������ļ���
*@return     void
*/void Reshaper::SaveBinVerts(const char * filename, const std::string & path, const std::vector<std::string>& files)
{
	Eigen::MatrixXd vertices;
	vertices.resize(VERTS * 3, files.size());
	//cout << "vertices.shape = " << vertices.rows() << " " << vertices.cols() << endl;
	int k = 0;
	for (auto file : files)
	{
		Eigen::Matrix3Xd vv;
		Eigen::Matrix3Xi ff;
		meshio::ReadObj(path + file, vv, ff);
		int cnt = 0;
		//std::cout << vv.cols() << " " << vv.rows() << std::endl;
		// be careful with the cols and rows
		for (int i = 0; i < vv.cols(); ++i)
		{
			for (int j = 0; j < vv.rows(); ++j)
			{
				vertices(cnt++, k) = vv(j, i);
			}
		}
		//cout << "cnt = " << file << endl;
		k++;
	}
	binaryio::WriteMatrixBinaryToFile(filename, vertices);
	//cout << "ok" << endl;
}

/*!
*@brief  �Զ����Ʊ�����Ƭ��Ϣ
*@param[out]
*@param[in]  const char * filename  ������ļ���
*@param[in]  const std::string & path  ģ�����ݼ�·��
*@param[in]  const std::vector<std::string> & files  ģ���ļ���
*@return     void
*/void Reshaper::SaveBinFaces(const char * filename, const std::string & path, const std::vector<std::string>& files)
{
	Eigen::Matrix3Xd vv;
	Eigen::Matrix3Xi ff;
	meshio::ReadObj(path + files[0], vv, ff);
	binaryio::WriteMatrixBinaryToFile(filename, ff);
	std::cout << "Save facets done!" << std::endl;
}

/*!
*@brief  �Զ����Ʊ�������ģ�͵Ķ�����Ϣ����Ƭ��Ϣ
*@param[out]
*@return     void
*/void Reshaper::SaveVertFacetInBin()
{
	std::string trainModelPath = DATASET_PATH;
	std::vector<std::string> trainFiles = GetFiles(trainModelPath + "*");

	SaveBinVerts((BIN_DATA_PATH + "vertex").c_str(), trainModelPath, trainFiles);

	// F (3, 12500)
	SaveBinFaces((BIN_DATA_PATH + "facets").c_str(), trainModelPath, trainFiles);

	Eigen::MatrixXd verts;
	Eigen::Matrix3Xi facets;

	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "vertex").c_str(), verts);
	binaryio::ReadMatrixBinaryFromFile((BIN_DATA_PATH + "facets").c_str(), facets);

	std::cout << facets.cols() << std::endl;
	std::cout << "Save facets verts done!" << std::endl;
}


