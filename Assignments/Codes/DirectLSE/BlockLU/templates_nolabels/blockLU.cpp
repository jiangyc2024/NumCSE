//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <cstdlib>
#include <Eigen/Dense>

using namespace Eigen;

/* @brief Solve the system Ry=c
 * for the upper triangular matrix R
 * This could help you in your implementation
 * of solve_LSE()
 * \param[in] R nxn regular, upper triangular matrix
 * \param[in] c n dim vector
 * \return y n dim result vector
 */
VectorXd solve_R(const MatrixXd& R, const VectorXd& c)
{
    int n = R.rows();
    assert(n == R.cols() && n == c.size() &&
           "Input dimensions must agree");
    // Initialize
    VectorXd y(n);
    // Implementing this function could help you in solve_LSE()
    return y;
}

/* @brief Solve the System Ax=b
 * for A << R,              v,
 *          u.transpose(),  0;
 * \param[in] R nxn regular, upper triangular matrix
 * \param[in] v n dim vector
 * \param[in] u n dim vector
 * \param[in] b n+1 dim vector
 * \return x n+1 dim result vector
 */
VectorXd solve_LSE(const MatrixXd& R,
               const VectorXd& v,
               const VectorXd& u,
               const VectorXd& b)
{
    unsigned n = R.rows();
    assert(R.cols() == n && "R has to be square");
    assert(n == v.size() && n == u.size() && n+1 == b.size()
           && "Input dimensions must agree");
    // Initialize
    VectorXd y(n+1), x(n+1);
    // Solve the LSE using LU-decomposition and the expression
    // for L and U that you derived
    return x;
}

int main()
{    
    // Vectors for testing
    unsigned n = 10;
    VectorXd v,u,b;
    u = v = VectorXd::Random(n);
    b = VectorXd::Random(n+1);
    // Upper triangular matrix
    MatrixXd R(n,n);
    for (unsigned i = 0; i < n; ++i)
    {
        for (unsigned j = i; j < n; ++j)
        {
            R(i,j) = rand(); //Bad RNG, but sufficient here
        }
    }
    R /= RAND_MAX;  //"norm" R for numerical stability 
    // Build matrix A for Eigensolver
    MatrixXd A(n+1,n+1);
    A << R,            v,
        u.transpose(), 0;
    
    double error = (solve_LSE(R,v,u,b) - A.colPivHouseholderQr().solve(b)).norm();
    if (error > 1e-8)
    {
        std::cout << "solve_LSE() returns a different solution than Eigen" << std::endl;
    } else
    {
        std::cout << "solve_LSE() and Eigen get the same result" << std::endl;
    }
}
