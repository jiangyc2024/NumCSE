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
/* SAM_LISTING_BEGIN_2 */
VectorXd solve_R(const MatrixXd& R, const VectorXd& c)
{
    int n = R.rows();
    assert(n == R.cols() && n == c.size() &&
           "Input dimensions must agree");
    // Initialize
    VectorXd y(n);
#if SOLUTION    
    // Since R is upper triangular, we can solve by backwards substitution
    for (int i = n-1; i >= 0; --i)
    {
        y(i) = c(i);
        for (int j = n-1; j > i ; --j)
        {
            y(i) -= R(i,j) * y(j);
        }
        y(i) /= R(i,i);
    }
#else
    // Implementing this function could help you in solve_LSE()
#endif
    return y;
}
/* SAM_LISTING_END_2 */

/* @brief Solve the System Ax=b
 * for A << R,              v,
 *          u.transpose(),  0;
 * \param[in] R nxn regular, upper triangular matrix
 * \param[in] v n dim vector
 * \param[in] u n dim vector
 * \param[in] b n+1 dim vector
 * \return x n+1 dim result vector
 */
/* SAM_LISTING_BEGIN_1 */
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
#if SOLUTION    
    // Solve the system Ax=b by LU-Decomposition.
    // Solve Ly = b through forward substitution.
    // Due to the special structure of our L,
    // the first n entries of y are easy:
    y.head(n) = b.head(n);
    // The last element of y is given by $y_n = b_n - u^T\mathbf{R}^{-1}y_{0...n-1}$
    y(n) = b(n) - u.transpose() * solve_R(R, y.head(n));
    
    // Solve Ux = y by backward substitution.
    // First we build U
    MatrixXd U(n+1,n+1);
    U << R,                             v,
         VectorXd::Zero(n).transpose(), -u.transpose()*solve_R(R,v);
    // Then we solve Ux = y
    x = solve_R(U,y);
    // \iffalse Latex comment-delimiter so this part doesn't appear
    // in the solution and screws up formatting
    // Note this could be done using less memory by not constructing
    // U, doing the first step of the back substitution "by hand" and then
    // calling solve_R()
    // VectorXd x(n+1);
    // x(n) = y(n) / (-u.transpose()*solve_R(R,v));
    // x.head(n) = solve_R(R,y.head(n) - v*x(n));
    //\fi
#else
    // Solve the LSE using LU-decomposition and the expression
    // for L and U that you derived
#endif
    return x;
}
/* SAM_LISTING_END_1 */

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
