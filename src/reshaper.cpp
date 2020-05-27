#include "reshaper.h"


Reshaper::Reshaper(){}

Reshaper::~Reshaper(){}


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
	for (auto file : files) {
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


