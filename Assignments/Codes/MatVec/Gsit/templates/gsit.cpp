#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

using namespace Eigen;

/* \brief Use symmetric Gauss-Seidel iterations to solve the system Ax = b
 * \param[in] A system matrix to be decompsed (L + D + U)
 * \param[in] b r.h.s. vector
 * \param[in,out] x initial guess and last iterate (approximated solution)
 * \param[in] rtol relative tolerance for termination criteria
 */
void GSIt(const MatrixXd & A, const VectorXd & b,
          VectorXd & x, double rtol) {
    /* SAM_LISTING_BEGIN_1 */
    // TODO: Implement Gauss-Seidel iteration
    /* SAM_LISTING_END_1 */

    return;
}

int main(int, char**) {
    unsigned int n = 9;

    MatrixXd A = MatrixXd::Zero(n,n);
    for(unsigned int i = 0; i < n; ++i) {
        if(i > 0) A(i,i-1) = 2;
        A(i,i) = 3;
        if(i < n-1) A(i,i+1) = 1;
    }
    VectorXd b = VectorXd::Constant(n, 1);

    //// Subproblem b: test GSIt
    std::cout << "--> Test GSIt" << std::endl;

    /* SAM_LISTING_BEGIN_2 */
    // TODO: test the code GSIt using the given data
    /* SAM_LISTING_END_2 */
}
