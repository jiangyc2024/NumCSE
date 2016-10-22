#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

/*
 * 2 dimensional sinus transform on matrix Y
 */
void sinetransform2d(const Eigen::MatrixXd& Y, Eigen::MatrixXd& S)
{
	int m = Y.rows();
	int n = Y.cols();

	Eigen::VectorXcd c;
	Eigen::FFT<double> fft;

	Eigen::MatrixXcd C(2*m+2,n);
	C.row(0) = Eigen::VectorXcd::Zero(n);
	C.middleRows(1, m) = Y.cast<std::complex<double>>();
	C.row(m+1) = Eigen::VectorXcd::Zero(n);
	C.middleRows(m+2, m) = -Y.colwise().reverse().cast<std::complex<double>>();

	// fft on each column of C
	for (int i=0; i<n; ++i)
	{
		fft.fwd(c,C.col(i));
		C.col(i) = c;
	}

	std::cout << C << std::endl << std::endl;

	std::complex<double> i(0,1);
	C.middleRows(1,m) = i*C.middleRows(1,m)/2.;

	Eigen::MatrixXcd C2(2*n+2,m);
	C2.row(0) = Eigen::VectorXcd::Zero(m);
	C2.middleRows(1,n) = C.middleRows(1,m).transpose();
	C2.row(n+1) = Eigen::VectorXcd::Zero(m);
	C2.middleRows(n+2, n) = -C.middleRows(1,m).transpose().colwise().reverse();

	// fft on each column of C2
	for (int i=0; i<m; ++i)
	{
		fft.fwd(c,C2.col(i));
		C2.col(i) = c;
	}

	S = (i*C2.middleRows(1,n).transpose() / 2.).real();
}

/*
void sinft2d(Y)
{
	[m,n] = size(Y);
	C = fft([zeros(1,n); Y;...
			 zeros(1,n);...
			 -Y(end:-1:1,:)]);
	C = i*C(2:m+1,:)'/2;
	C = fft([zeros(1,m); C;...
			 zeros(1,m);...
			 -C(end:-1:1,:)]);
	C= i*C(2:n+1,:)'/2;

	return C;
}*/
