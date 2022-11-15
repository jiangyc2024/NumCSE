#include <Eigen/Sparse>
#include <fstream>
#include <iostream>
#include <string>

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
int loadGraphMarketMatrix(Eigen::SparseMatrix<type> &A, const std::string & path)
{
	std::ifstream matFile;
	matFile.open(path.c_str());
	int code = 0;

	if (matFile.is_open())
	{
		std::string line;
		std::vector<Eigen::Triplet<type>> triplets;

		bool read_sizes= false;
		Eigen::Index M = 0;
		Eigen::Index N = 0;
		Eigen::Index count = 0;

		while (std::getline(matFile, line))
		{
			if (line[0] == '%') {
				
				continue;
			}

			std::stringstream strstream(line);

			if (! read_sizes)
			{
				strstream >> M >> N >> count;

				if (M > 0 && N > 0 && count > 0)
				{
					read_sizes = true;
					triplets.reserve(count);
				}
			}	
			else
			{
				Eigen::Index i = -1;
				Eigen::Index j = -1;
				strstream >> i >> j;
				if (i >= 0 && j >= 0 && i < M && j < N)
				{
					triplets.emplace_back(i,j,1);
				}
			}

		}

		A.resize(M, N);
		A.setFromTriplets(triplets.begin(), triplets.end());

		code = 1;
	}
	else {

		std::cerr << "Error reading file " << path << std::endl;
	}

	matFile.close();

	std::quick_exit(EXIT_FAILURE);
	return code;
}


template <class type>
int loadGraphMarketMatrix(Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic> &A, std::string path)
{
	
	int code = 0;
	Eigen::SparseMatrix<type> GSparse;
 	if (loadGraphMarketMatrix(GSparse, path))
	{	
		A = Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic>(GSparse);
		code = 1;
	}
	
	return code;
}
