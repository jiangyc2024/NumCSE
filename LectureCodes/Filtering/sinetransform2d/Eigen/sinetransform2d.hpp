#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

/* SAM_LISTING_BEGIN_0 */
void sinetransform2d(const Eigen::MatrixXd& Y, Eigen::MatrixXd& S)
{
	int m = Y.rows();
	int n = Y.cols();

	Eigen::VectorXcd c;
	Eigen::FFT<double> fft;
	std::complex<double> i(0,1);

	Eigen::MatrixXcd C(2*m+2,n);
	C.row(0) = Eigen::VectorXcd::Zero(n);
	C.middleRows(1, m) = Y.cast<std::complex<double>>();
	C.row(m+1) = Eigen::VectorXcd::Zero(n);
	C.middleRows(m+2, m) = -Y.colwise().reverse().cast<std::complex<double>>();

	// FFT on each column of C -- Eigen::fft only operates on vectors
	for (int i=0; i<n; ++i)
	{
		fft.fwd(c,C.col(i));
		C.col(i) = c;
	}
	
	C.middleRows(1,m) = i*C.middleRows(1,m)/2.;

	Eigen::MatrixXcd C2(2*n+2,m);
	C2.row(0) = Eigen::VectorXcd::Zero(m);
	C2.middleRows(1,n) = C.middleRows(1,m).transpose();
	C2.row(n+1) = Eigen::VectorXcd::Zero(m);
	C2.middleRows(n+2, n) = -C.middleRows(1,m).transpose().colwise().reverse();

	// FFT on each column of C2 -- Eigen::fft only operates on vectors
	for (int i=0; i<m; ++i)
	{
		fft.fwd(c,C2.col(i));
		C2.col(i) = c;
	}

	S = (i*C2.middleRows(1,n).transpose() / 2.).real();
}
/* SAM_LISTING_END_0 */
