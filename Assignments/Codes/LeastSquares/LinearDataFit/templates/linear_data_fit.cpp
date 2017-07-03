#include <iostream>
#include <cmath>

#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace Eigen;

/* SAM_LISTING_BEGIN_1 */
MatrixXd make_A(const VectorXd &b) {
	size_t n = b.size();
	MatrixXd A(n, 4);
	// TODO: construct fitting matrix A
	return A;
}
/* SAM_LISTING_END_1 */


/* SAM_LISTING_BEGIN_2 */
VectorXd data_fit_normal(const MatrixXd &A, const VectorXd &b) {
	// TODO: solve least-squares problem using multiplication with transpose
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
VectorXd data_fit_qr(const MatrixXd &A, const VectorXd &b) {
	// TODO: solve least-squares problem using QR-Decomposition
}
/* SAM_LISTING_END_3 */


/* SAM_LISTING_BEGIN_4 */
int main() {
	// TODO: create a lin-log plot using the mgl::Figure class

	// TODO: analyze errors
}
/* SAM_LISTING_END_4 */
