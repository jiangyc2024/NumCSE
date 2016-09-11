# pragma once
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
using Eigen::VectorXcd;

VectorXcd pconvfft(const VectorXcd& u, const VectorXcd& x) {
  Eigen::FFT<double> fft;
  VectorXcd tmp = ( fft.fwd(u) ).cwiseProduct( fft.fwd(x) );
  return fft.inv(tmp);
}
