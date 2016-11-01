#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

using namespace Eigen;

/* @brief Compute the matrix $C$ from $A$
 * @param[in] A An $n \times n$ matrix
 * @param[out] C The $(n^2) \times (n^2)$ matrix
 * from $A\otimes I+I\otimes A$
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd polyMult_naive(const VectorXd & u, const VectorXd & v)
{
    // Initialization
    unsigned m = u.size() - 1; // Degree of polynomial
    unsigned n = v.size() - 1; // Degree of polynomial
    unsigned dim = max(m, n);
    
    VectorXd u_tmp = u;
    u_tmp.conservativeResize(dim);
    VectorXd v_tmp = v;
    v_tmp.conservativeResize(dim);
    
    VectorXd uv(m+n+1);

#if SOLUTION
    for(unsigned i=0; i<uv.size(); ++i) {
		unsigned fst = max(0, i - n);
		unsigned lst = min(i, n);
		for(unsigned j=fst; j<=lst; ++j) {
			uv(i) += u_tmp(j) * v_tmp(i-j);
		}
	}
#else // TEMPLATE
    // TODO: compute $C$ from $A$
#endif // TEMPLATE

	return uv;
}
/* SAM_LISTING_END_0 */

/* @brief Solve the Lyapunov system
 * @param[in] A An $n \times n$ matrix
 * @param[out] X The $n \times n$ solution matrix
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd polyMult_fast(const VectorXd & u, const VectorXd & v)
{
    // Initialization
    unsigned m = u.size() - 1; // Degree of polynomial
    unsigned n = v.size() - 1; // Degree of polynomial
    unsigned dim = max(m, n);
    
    VectorXd u_tmp = u;
    u_tmp.conservativeResize(u.size() + n);
    VectorXd v_tmp = v;
    v_tmp.conservativeResize(v.size() + m);

#if SOLUTION
	Eigen::FFT<double> fft;
	VectorXcd tmp = ( fft.fwd(u_tmp) ).cwiseProduct( fft.fwd(v_tmp) );
	VectorXd uv = fft.inv(tmp).real();
#else // TEMPLATE
    // TODO: compute $C$ from $A$
#endif // TEMPLATE

	return uv;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
int main() {
    // Initialization
    unsigned int n = 5;
    MatrixXd A(n,n);
    A << 10, 2, 3, 4, 5,
         6, 20, 8, 9, 1,
         1, 2, 30, 4, 5,
         6, 7, 8, 20, 0,
         1, 2, 3, 4, 10;

    // Test 'buildC'
    SparseMatrix <double> C = buildC(A);
    std::cout << "C = " << C << std::endl;

    // Test 'solveLyapunov'
    MatrixXd X(n,n);
    solveLyapunov(A,X);
    std::cout << "X = " << X << std::endl;

    // Verify the solution if you obtain zero
    MatrixXd I = MatrixXd::Identity(n,n);
    std::cout << "Correct if close to 0: "
              << (A*X + X*A.transpose() - I).norm()
              << std::endl;
}
/* SAM_LISTING_END_2 */
