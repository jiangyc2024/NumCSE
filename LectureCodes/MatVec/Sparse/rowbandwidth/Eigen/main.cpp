///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "bandwidth.hpp"

using namespace Eigen;

int main () {
	int n = 7;
	typedef Triplet<double> triplet_t;
	std::vector<triplet_t> triplets;
	
	// example in lecture document ex:env
	for(int i = 0; i < n; ++i)
		triplets.push_back( {i,i,1.0} );//diagonal
	triplets.push_back( {0,2,1.0} );
	triplets.push_back( {1,4,1.0} );
	triplets.push_back( {2,0,1.0} );
	triplets.push_back( {2,6,1.0} );
	triplets.push_back( {3,4,1.0} );
	triplets.push_back( {3,6,1.0} );
	triplets.push_back( {4,1,1.0} );
	triplets.push_back( {4,3,1.0} );
	triplets.push_back( {4,5,1.0} );
	triplets.push_back( {5,4,1.0} );
	triplets.push_back( {6,2,1.0} );
	triplets.push_back( {6,3,1.0} );
	SparseMatrix<double> A(n,n);
	A.setFromTriplets(triplets.begin(), triplets.end());
	A.makeCompressed();
	
	MatrixXd B(A);
	
	std::cout << A << std::endl;
	std::cout << B << std::endl;
	
	std::cout << "Rowbandwidth\n" << rowbandwidth(A) << std::endl << std::endl;
	std::cout << "Rowbandwidth\n" << rowbandwidth(B) << std::endl << std::endl;
	
	std::cout << "colbandwidth\n" << colbandwidth(A) << std::endl << std::endl;
	std::cout << "colbandwidth\n" << colbandwidth(B) << std::endl << std::endl;
	return 0;
}

