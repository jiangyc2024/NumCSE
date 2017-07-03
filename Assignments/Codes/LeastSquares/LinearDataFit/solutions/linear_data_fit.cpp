#include <iostream>
#include <cmath>

#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace Eigen;

/* SAM_LISTING_BEGIN_1 */
MatrixXd make_A(const VectorXd &b) {
	size_t n = b.size();
	MatrixXd A(n, 4);
	for (size_t i = 0; i < n; i++) {
		double t = b(i);
		A(i, 0) = 1.0 / t;
		A(i, 1) = 1.0 / (t*t);
		A(i, 2) = std::exp(-(t-1));
		A(i, 3) = std::exp(-2*(t-1));
	}
	return A;
}
/* SAM_LISTING_END_1 */


/* SAM_LISTING_BEGIN_2 */
VectorXd data_fit_normal(const MatrixXd &A, const VectorXd &b) {
	auto At = A.transpose();
	auto AtA = At * A;
	auto Atb = At * b;
	return AtA.ldlt().solve(Atb);
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
VectorXd data_fit_qr(const MatrixXd &A, const VectorXd &b) {
	return A.colPivHouseholderQr().solve(b);
}
/* SAM_LISTING_END_3 */


/* SAM_LISTING_BEGIN_4 */
int main() {
	auto t = VectorXd::LinSpaced(10, 0.1, 1.0);
	auto A = make_A(t);

	VectorXd f(10);
	f << 100. , 34. , 17. , 12. , 9. , 6. , 5. , 4. , 4. , 2.;

	auto gamma1 = data_fit_normal(A, f);
	auto gamma2 = data_fit_qr(A, f);

	auto y1 = A*gamma1;
	auto y2 = A*gamma2;

	auto tl = VectorXd::LinSpaced(91, 0.1, 1.0);
	auto Al = make_A(tl);
	auto yl1 =  Al*gamma1;
	auto yl2 =  Al*gamma2;


	mgl::Figure fig1;
	fig1.setlog(false, true);
	fig1.plot(tl, yl1, "r").label("normal equation");
	fig1.plot(tl, yl2, "b").label("qr fitting");
	fig1.plot(t, f, "k*").label("data set");
	fig1.xlabel("t");
	fig1.ylabel("y");
	fig1.legend(1, 1);
	fig1.save("fitted.eps");


	VectorXd err1 = (y1-f);
	err1 = err1.cwiseProduct(err1);
	VectorXd err2 = (y2-f);
	err2 = err2.cwiseProduct(err2);

	mgl::Figure fig2;
	fig2.setlog(false, true);
	fig2.plot(t, err1, "r*").label("normal equation");
	fig2.plot(t, err2, "bo").label("qr fitting");
	fig2.xlabel("t");
	fig2.ylabel("fitting error");
	fig2.legend(1, 1);
	fig2.save("errors.eps");

	std::cout << (gamma1 - gamma2) << std::endl;
	std::cout << "L2-Norms: " << std::sqrt(err1.sum()) << " " << std::sqrt(err2.sum()) << std::endl;
	std::cout << "Difference in L2-Norms: " << (std::sqrt(err1.sum()) - std::sqrt(err2.sum())) << std::endl;

	auto cond = [](MatrixXd A){
		JacobiSVD<MatrixXd> svd(A);
		auto sigma = svd.singularValues();
		return sigma[0]/sigma[sigma.size()-1];
	};

	std::cout << "cond(A)" << cond(A) << std::endl;
	std::cout << "cond(A^t A)" << cond(A.transpose() * A) << std::endl;

}
/* SAM_LISTING_END_4 */
