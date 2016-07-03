#include <string>
#include <iostream>
#include <fstream>

/*
 * loads a matrix in market format that represents a graph (only indices but no values given.
 * For example:
 * 10 10 2
 * 0 5
 * 4 3
 *
 * This is unfortunately not supported by the Eigen::loadMarket function in unsupported/Eigen/SparseExtra>
 */

template <class type>
int loadGraphMarketMatrix(Eigen::SparseMatrix<type> &A, std::string path)
{
	std::ifstream matFile;
	matFile.open(path.c_str());


	if (matFile.is_open())
	{
		std::string line;
		std::vector<Eigen::Triplet<type>> triplets;

		int readsizes = 0;
		int M, N, count;

		while (std::getline(matFile, line))
		{
			if (line[0] == '%') continue;
			std::stringstream strstream(line);

			if (readsizes == 0)
			{
				strstream >> M >> N >> count;

				if (M > 0 && N > 0 && count > 0)
				{
					readsizes = 1;
					triplets.reserve(count);
				}
			}	
			else
			{
				int i(-1), j(-1);
				strstream >> i >> j;
				if (i >= 0 && j >= 0 && i < M && j < N)
				{
					triplets.push_back(Eigen::Triplet<type>(i,j,1));
				}
			}

		}

		A.resize(M, N);
		A.setFromTriplets(triplets.begin(), triplets.end());

		return 1;
	}
	else std::cerr << "Error reading file " << path << std::endl;

	matFile.close();

	exit(EXIT_FAILURE);
	return 0;
}


template <class type>
int loadGraphMarketMatrix(Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic> &A, std::string path)
{
	Eigen::SparseMatrix<type> GSparse;
 	if (loadGraphMarketMatrix(GSparse, path))
	{	
		A = Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic>(GSparse);
		return 1;
	}
	else return 0;
}
