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

using namespace Eigen;

/* \brief Solve system Ax = b with optimal complexity O(n)
 * \param[in] c Entries for matrix $A$
 * \param[in] b r.h.s. vector
 * \param[in] i0 index
 * \param[in] j0 index
 * \return Solution x, s.t. Ax = b
 */
VectorXd solveLSE(const VectorXd & c, const VectorXd & b,
                  int i0, int j0) {
    assert(c.size() == b.size()-1
           && "Size mismatch!");
    assert(i0 > j0 && "i0 must be bigger than j0!");

    VectorXd ret(b.size());

    ret(b.size()-1) = b(b.size()-1);

    for(int i = b.size()-2; i >= 0; --i) {
        if(i == i0) {
            double piv = 1;
            double rhs = b(i);
            int sign = -1;
            for(int j = j0; j < i0; ++j) {
                rhs -= piv * b(j);
                piv *= sign*c(j);
                sign *= -1;
                std::cout << "a" << c(j);
            }
            std::cout << "a" << piv;
            ret(i) = (rhs - ret(i+1)*c(i)) / (1.-piv);
        } else {
            ret(i) = b(i) - ret(i+1)*c(i);
        }
    }

    return ret;
}

int main(int, char**) {
    // Size of matrix
    unsigned int n = 15;

    // Vector of triplets
    std::vector< Triplet<double> > triplets;

    unsigned int ntriplets = 2*n;

    // Two random indices
    unsigned int i0 = 6, j0 = 4;

    //// PROBLEM a
    std::cout << "--> PROBLEM a:" << std::endl;

    // Reserve space for triplets
    triplets.reserve(ntriplets);

    // Used for SolveLSE
    VectorXd c(n-1);

    // Build triplets vector
    for(unsigned int i = 0; i < n; ++i) {
        // Used for SolveLSE
        if(i < n-1) {
            c(i) = i;
        }

        // Build matrix using triplet format
        triplets.push_back(
                           Triplet<double>(i, i, 1)
                           );
        if(i < n-1) {
            triplets.push_back(
                               Triplet<double>(i, i+1, c(i))
                               );
        }
    }
    triplets.push_back(
                       Triplet<double>(i0, j0, 1)
                       );

    // Build system
    SparseMatrix<double> A(n,n);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();

    // Solve system
    SparseLU< SparseMatrix<double> > splu;
    splu.analyzePattern(A);
    splu.factorize(A);

    // Compute and output error
    VectorXd b = VectorXd::Random(n);
    std::cout << "Error:" << std::endl
              << (splu.solve(b) - solveLSE(c, b, i0, j0)).norm()
              << std::endl;
}
