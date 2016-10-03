//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iostream>

#include <vector>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using Triplet = Eigen::Triplet<double>;
using Triplets = std::vector<Triplet>;

using Vector = Eigen::VectorXd;
using Matrix = Eigen::SparseMatrix<double>;

//! \brief Solve system Ax = b with optimal complexity O(n)
//! \param[in] c Entries for matrix A
//! \param[in] b r.h.s. vector
//! \param[in] i0 index
//! \param[in] j0 index
//! \return Solution x, s.t. Ax = b
Vector solveLSE(const Vector & c, const Vector & b, unsigned int i0, unsigned int j0) {
    assert(c.size() == b.size()-1 && "Size mismatch!");
    
    Vector ret(b.size());
    ret(b.size()-1) = b(b.size()-1);
    for(int i = b.size()-2; i >= 0; --i) {
        if(i == (int) i0) {
            double piv = 1;
            double rhs = b(i);
            int sign = -1;
            for(unsigned int j = j0; j < i0; ++j) {
                rhs -= piv * b(j);
                piv *= sign*c(j);
                sign *= -1;
                std::cout << "a"<< c(j);
            }
            std::cout << "a"<< piv;
            ret(i) = (rhs - ret(i+1)*c(i)) / (1.-piv);
        } else {
            ret(i) = b(i) - ret(i+1)*c(i);
        }
    }
    
    return ret;
}

int main(int, char**) {
    unsigned int n = 15;
    
    Triplets triplets;
    
    unsigned int ntriplets = 2*n;
    
    unsigned int i0 = 6, j0 = 4;
    
    //// PROBLEM 6a
    std::cout << "*** PROBLEM 6a:" << std::endl;
    
    // reserve space
    triplets.reserve(ntriplets);
    
    Vector c(n-1);
    
    // Build triplets vecotr
    for(unsigned int i = 0; i < n; ++i) {
        if(i < n-1) c(i) = i;
        triplets.push_back(Triplet(i,i,1));
        if(i < n-1) triplets.push_back(Triplet(i,i+1,c[i]));
    }
    triplets.push_back(Triplet(i0,j0,1));
    
    // Solve system
    Matrix A(n,n);
    A.setFromTriplets(triplets.begin(), triplets.end());
    
    A.makeCompressed();
    Eigen::SparseLU<Matrix> splu;
    splu.analyzePattern(A); 
    splu.factorize(A);
    
    Vector b = Vector::Random(n);
    std::cout << "Eigen Solve:" << std::endl << splu.solve(b) << std::endl;
    std::cout << "My Solve:" << std::endl << solveLSE(c,b,i0,j0) << std::endl;
}
