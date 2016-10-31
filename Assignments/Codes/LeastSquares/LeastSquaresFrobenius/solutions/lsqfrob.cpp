#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse> // FIX bug in Eigen
#include <unsupported/Eigen/KroneckerProduct>

using namespace Eigen;

//#define DECOMP partialPivLU
//#define DECOMP fullPivLU
//#define DECOMP llt
#define DECOMP ldlt

/* @brief Solve the minimization problem $\argmin |M|_{F}$
 * @param[in] z An $n$-dimensional vector constraining the solution $M^*$ via $M^*z = g$
 * @param[in] g An $n$-dimensional vector constraining the solution $M^*$ via $M^*z = g$
 * @param[out] M^* The solution to the minimization problem
 */
/* SAM_LISTING_BEGIN_1 */
MatrixXd min_frob(const VectorXd & z, const VectorXd & g) {
    assert(z.size() == g.size() && "Size mismatch!");

    unsigned int n = g.size();
    
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

    MatrixXd Mstar = min_frob(z, g); // $M^*$
    MatrixXd M = g*z.transpose() / z.squaredNorm(); // $M$

    std::cout << "Norm: "
              << (Mstar - M).norm()
              << std::endl;
/* SAM_LISTING_END_2 */
}
