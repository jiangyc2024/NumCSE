#include <iostream>

#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

using namespace Eigen;

/* SAM_LISTING_BEGIN_1 */
MatrixXd min_frob(const VectorXd & z, const VectorXd & g) {
    assert(z.size() == g.size() && "Size mismatch!");

    unsigned int n = g.size();

    // Build temporary matrix $\mathbf{C}$.
    MatrixXd C = kroneckerProduct(MatrixXd::Identity(n,n),
                                  z.transpose());

    // Build system matrix and r.h.s.
    MatrixXd lhs(n*n+n,n*n+n);
    S << MatrixXd::Identity(n*n, n*n), C.transpose(),
         C                           , MatrixXd::Zero(n,n);
    VectorXd h(n*n+n);
    h << VectorXd::Zero(n*n), g;

    // Solve augmented system and return only head of solution
    // as a $n \times n$ matrix (discard the rest).
    return Map<const MatrixXd>(S.lu().solve(h).eval().data(),
                               n, n);
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

    MatrixXd Mstar = min_frob(z, g).transpose(); // $\mathbf{M}^*$
    MatrixXd M = g*z.transpose() / z.squaredNorm(); // $\mathbf{M}$

    std::cout << "Norm: "
              << (Mstar - M).norm()
              << std::endl;
    /* SAM_LISTING_END_2 */
}
