#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
void applyHouseholder(VectorXd &x, const VectorXd &v) {
    // TODO
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
template <typename Scalar>
void applyHouseholder(VectorXd &x,
                      const MatrixBase<Scalar> &V) {
    int n = V.rows();

    // TODO
}
/* SAM_LISTING_END_1 */

int main() {

}
