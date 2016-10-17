#include <iostream>
#include <cmath>

#include <Eigen/Dense>

using namespace Eigen;


/* SAM_LISTING_BEGIN_1 */
/*!
 * \brief circl_alg_fit
 * \param[in] x
 * \param[in] y
 * \return
 */
Vector3d circl_alg_fit(const VectorXd &x,
                       const VectorXd & y) {
    assert(x.size() == y.size() && "Size mismatch!");

    unsigned int n = x.size();

    // TODO: find center/radius using algebriac method
}
/* SAM_LISTING_END_1 */



/* SAM_LISTING_BEGIN_2 */
Vector3d circl_geo_fit(const VectorXd &x, const VectorXd & y) {

}
/* SAM_LISTING_END_2 */



Vector3d circl_geo_fit(const VectorXd &x, const VectorXd & y) {
//    A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b)
//            A.colPivHouseholderQr().solve(b) << endl;
//    A.transpose() * A).ldlt().solve(A.transpose() * b) << endl;
}

int main(int argc, char **argv) {

}
