//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

#include "ellpack.hpp"

using namespace Eigen;


 /* SAM_LISTING_BEGIN_5 */
int main(int, char**) {
    // Vector of triplets
    Triplets triplets;

    // Data
    unsigned int m = 3, n = 6;
    unsigned int ntriplets = 9;

    // Reserve space for triplets
    triplets.reserve(ntriplets);

    // Fill in some triplets
    triplets.push_back(Triplet_new(1,2,4));
    triplets.push_back(Triplet_new(0,0,5));
    triplets.push_back(Triplet_new(1,2,6));  //repeated index
    triplets.push_back(Triplet_new(2,5,7));
    triplets.push_back(Triplet_new(0,4,8));
    triplets.push_back(Triplet_new(0,0,1));  //repeated index
    triplets.push_back(Triplet_new(1,3,9));
    triplets.push_back(Triplet_new(2,2,10));
    triplets.push_back(Triplet_new(1,3,2));  //repeated index
    triplets.push_back(Triplet_new(2,1,11));
    triplets.push_back(Triplet_new(1,0,12));

    // Build Ellpack matrix
    EllpackMat E(triplets, m, n);

    //// TEST
    std::cout << " ------------- Test: print  E ---------------- " << std::endl;
    for (int i = 0; i < m; ++i ) {
      for (int j = 0; j < n; ++j ) {
        std::cout << std::setw(4) << E(i,j) << " ";
      }
      std::cout << std::endl;
    }


    std::cout << " ------------- Test of y = E*x --------------- " << std::endl;
    Vector x(6);
    x << 4,5,6,7,8,9;
    
    Vector Ex = Vector::Zero(m);
    E.mvmult(x, Ex);
    std::cout << "Ellpack E*x =" << std::endl << Ex << std::endl;

    std::cout << " ------------- Test of y = A^t*x ------------- " << std::endl;
    x.resize(3);
    x << 1,2,3;
    
    Vector Etx = Vector::Zero(n);
    E.mtvmult(x, Etx);
    std::cout << "Ellpack E^t*x =" << std::endl << Etx << std::endl;
}
/* SAM_LISTING_END_5 */
