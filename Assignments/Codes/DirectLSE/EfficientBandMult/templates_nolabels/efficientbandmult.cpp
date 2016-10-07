//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace Eigen;

/* @brief Compute $y = A*x$ with A banded matrix with diagonal structure
 * \param[in] a An $(n-1)$-dimensional vector for the upper diagonal
 * \param[in] b An $(n-2)$-dimensional vector for the second lower diagonal
 * \param[in] x An $n$-dimensional vector for $Ax = y$
 * \param[out] y The $n$-dimensional vector $y = Ax$
 */
template <class Vector>
void multAx(const Vector & a, const Vector & b, const Vector & x, Vector & y) {
    unsigned int n = x.size();
    if( a.size() < n - 1 || b.size() < n - 2 ) {
        // Size check, error if do not match
        std::cerr << "Error: size mismatch!" << std::endl;
        return;
    }
    y.resize(n);
    
    // TODO: compute $y = A*x$ with $A$ banded matrix
}

/* @brief Solve $r = A*x$ with $A$ banded matrix with upper triangular sparse structure
 * \param[in] a An $(n-1)$-dimensional vector for the upper diagonal
 * \param[in] r An $n$-dimensional vector for $Ax = r$
 * \param[out] x The $n$-dimensional vector from $Ax = r$
 */
template <class Vector>
void solvelseAupper(const Vector & a, const Vector & r, Vector & x) {
    // Set up dimensions
    unsigned int n = r.size();
    x.resize(n);
    
    // TODO: solve system $r = A*x$ with $A$ upper triangular matrix
}

/* @brief Solve $r = A*x$ with $A$ banded matrix using Gaussian elimination (no pivot)
 * \param[in] a An $(n-1)$-dimensional vector for the upper diagonal
 * \param[in] b An $(n-2)$-dimensional vector for the second lower diagonal
 * \param[in] r An $n$-dimensional vector for $Ax = r$
 * \param[out] x The $n$-dimensional vector from $Ax = r$
 */
template <class Vector>
void solvelseA(const Vector & a, const Vector & b, const Vector & r, Vector & x) {
    // Set up dimensions
    int n = r.size();
    Vector c(n-1,0);
    Vector d(n,  2);
    Vector y;
    x = r;
    y = r;
      
    // TODO: solve system $r = A*x$ using Gaussian elimination
}

/* @brief Solve $r = A*x$ with $A$ banded matrix using Eigen::SparseLU
 * \param[in] a An $(n-1)$-dimensional vector for the upper diagonal
 * \param[in] b An $(n-2)$-dimensional vector for the second lower diagonal
 * \param[in] r An $n$-dimensional vector for $Ax = r$
 * \param[out] x The $n$-dimensional vector from $Ax = r$
 */
template <class Vector>
void solvelseAEigen(const Vector & a, const Vector & b, const Vector & r, Vector & x) {
    // Set up dimensions
    typedef typename Vector::Scalar Scalar;
    unsigned int n = r.size();
    
    // Fill in matrix:
    // We reserve 3 nonzero entries per row for Gaussian fill-in
    SparseMatrix<Scalar> A(n,n);
    A.reserve(3);
    for(unsigned int i = 0; i < n; ++i) {
        A.insert(i,i) = 2;
        if(i < n-1) A.insert(i,i+1) = a(i);
        if(i >= 2)  A.insert(i,i-2) = b(i-2);
    }
    A.makeCompressed();
    
    // TODO: solve system $r = A*x$ using Eigen::SparseLU
}

int main() {
    unsigned int n = 9;
    // Compute with all three solvers
    std::cout << "*** Check that the solvers are correct" << std::endl;
    VectorXd a = VectorXd::Random(n-1);
    VectorXd b = VectorXd::Zero(n-2); // All 0s for upper diagonal structure
    VectorXd y = VectorXd::Random(n);
    VectorXd x;
    
    std::cout << "Original: " << y << std::endl;
    
    // Compute $y = A*A^{-1}*y$
    solvelseAupper(a,y,x);
    multAx(a,b,x,y);
    
    // Should be the same as before
    std::cout << "Upper: " << y << std::endl;
    
    // Random $b$ for own and Eigen-based solver
    b = VectorXd::Random(n-2);

    // Compute $y = A*A^{-1}*y$
    solvelseA(a,b,y,x);
    multAx(a,b,x,y);
    
    // Should be the same as before
    std::cout << "Own: " << y << std::endl;
    
    // Compute $y = A*A^{-1}*y$
    solvelseAEigen(a,b,y,x);
    multAx(a,b,x,y);
    
    // Should be the same as before
    std::cout << "Eigen: " << y << std::endl;
}
