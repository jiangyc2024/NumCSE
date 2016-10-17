#include <iostream>

#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

using namespace Eigen;

//#define DECOMP partialPivLU
//#define DECOMP fullPivLU
//#define DECOMP llt
#define DECOMP ldlt

/* SAM_LISTING_BEGIN_1 */
MatrixXd min_frob(const VectorXd & z, const VectorXd & g) {
    assert(z.size() == g.size() && "Size mismatch!");

    unsigned int n = g.size();
#if SOLUTION

    // Build temporary matrix $\mathbf{C}$.
    MatrixXd C = kroneckerProduct(MatrixXd::Identity(n,n),
                                  z.transpose());

    // Build system matrix and r.h.s.
    MatrixXd S(n*n+n,n*n+n);
    S << MatrixXd::Identity(n*n, n*n), C.transpose(),
         C                           , MatrixXd::Zero(n,n);
    VectorXd h(n*n+n);
    h << VectorXd::Zero(n*n), g;

    // Solve augmented system and return only head of solution
    // as a $n \times n$ matrix (discard the rest).
    return Map<const MatrixXd>(S.DECOMP().solve(h).eval().data(),
                               n, n).transpose();
    // It maybe possible to solve the system more efficiently
    // eployting the special structure of the matrix.
#else // TEMPLATE
    // TODO: solve minimization problem, return the matrix $\mathbf{M}^*$
    return MatrixXd::Zero(n,n);
#endif // TEMPLATE
}
/* SAM_LISTING_END_1 */

int main(int argc, char **argv) {

    unsigned int n = 100;
    if(argc > 1) {
        n = std::stoi(argv[1]);
    }

    /* SAM_LISTING_BEGIN_2 */
    VectorXd z = VectorXd::Random(n),
             g = VectorXd::Random(n);

#if SOLUTION
    MatrixXd Mstar = min_frob(z, g); // $\mathbf{M}^*$
    MatrixXd M = g*z.transpose() / z.squaredNorm(); // $\mathbf{M}$
#else // TEMPLATE
    // TODO: Compute matrices $M$ and $M^*$
    MatrixXd Mstar = MatrixXd::Zero(n,n); // $\mathbf{M}^*$
    MatrixXd M = MatrixXd::Zero(n,n); // $\mathbf{M}$
#endif // TEMPLATE

    std::cout << "Norm: "
              << (Mstar - M).norm()
              << std::endl;
    /* SAM_LISTING_END_2 */
}
