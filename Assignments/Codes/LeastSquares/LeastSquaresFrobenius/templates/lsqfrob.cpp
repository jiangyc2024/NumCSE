#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse> // FIX bug in Eigen
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
    // TODO: solve minimization problem, return the matrix $\mathbf{M}^*$
    return MatrixXd::Zero(n,n);
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

    // TODO: Compute matrices $M$ and $M^*$
    MatrixXd Mstar = MatrixXd::Zero(n,n); // $\mathbf{M}^*$
    MatrixXd M = MatrixXd::Zero(n,n); // $\mathbf{M}$

    std::cout << "Norm: "
              << (Mstar - M).norm()
              << std::endl;
    /* SAM_LISTING_END_2 */
}
