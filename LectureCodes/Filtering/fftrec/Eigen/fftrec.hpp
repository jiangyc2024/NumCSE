# include <iostream>
# include <complex>
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
using Eigen::VectorXcd;

void fftrec(const VectorXcd& y, VectorXcd& c) {
  const long n = y.size();
  if (n == 1) {  c = y; return; }
  if (n % 2 != 0) {
    std::cout << "size(y) must be even!\n";
    return;
  }

  Eigen::FFT<double> fft;
  VectorXcd c1( n/2 ), c2( n/2 ), yh( n/2 );
  for (long j = 0; j < n/2; ++j)  {
    yh(j) = y(2*j);
    c1 = fft.fwd(yh);
  }
  for (long l = 0; l < n/2; ++l) {
    yh(l) = y(2*l + 1);
    c2 = fft.fwd(yh);
  }
  std::complex<double> i(0,1); // imaginary unit
  c.resize(n);
  for (long k = 0; k < n; ++k) {
    c(k) = c1(k%(n/2)) + c2(k%(n/2))*std::exp(-2*M_PI/n*k*i);
  }
}
