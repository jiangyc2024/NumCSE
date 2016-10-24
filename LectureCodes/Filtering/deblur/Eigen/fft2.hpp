# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
using Eigen::VectorXcd;
using Eigen::MatrixXcd;

//! FFT for matrices
//  This just implements FFT for each column of X
MatrixXcd fft(const MatrixXcd& X) {
  const long m = X.rows(), n = X.cols();
  MatrixXcd Y(m, n);
  Eigen::FFT<double> fft;
  for (long j = 0; j < n; ++j) {
    VectorXcd Xj = X.col(j);
    Y.col(j) = fft.fwd(Xj);
  }
  return Y;
}

//! Inverse FFT for matrices 
//  This just implements inverse FFT for each column of X
MatrixXcd ifft(const MatrixXcd& X) {
  const long m = X.rows(), n = X.cols();
  MatrixXcd Y(m, n);
  Eigen::FFT<double> fft;
  for (long j = 0; j < n; ++j) {
    VectorXcd Xj = X.col(j);
    Y.col(j) = fft.inv(Xj);
  }
  return Y;
}

//! 2-dimensional FFT
//  Implementation based on: https://ch.mathworks.com/help/matlab/ref/fft2.html
MatrixXcd fft2(const MatrixXcd& X) {
  return fft(fft(X).transpose()).transpose();
}

//! 2-dimensional inverse FFT
MatrixXcd ifft2(const MatrixXcd& X) {
  return ifft(ifft(X).transpose()).transpose();
}
