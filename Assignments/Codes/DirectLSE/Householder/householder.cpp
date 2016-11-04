#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
void applyHouseholder(VectorXd &x, const VectorXd &v) {
#if SOLUTION
    double d = v.transpose()*x;
    x -= 2*v*d;
#else // TEMPLATE
    // TODO
#endif // TEMPLATE
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
template <typename Scalar>
void applyHouseholder(VectorXd &x,
                      const MatrixBase<Scalar> &V) {
    int n = V.rows();

#if SOLUTION
    for(unsigned int i = 0; i < n; ++i) {
        VectorXd v = V.col(i);
        v(i) = std::sqrt(1. - V.col(i).normSquared());

        x += 2.*v.dot(x) / (1. + 2.*v.dot(v)) * v;
    }
#else // TEMPLATE
    // TODO
#endif // TEMPLATE
}
/* SAM_LISTING_END_1 */

int main() {

}
