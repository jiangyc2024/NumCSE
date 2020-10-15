#ifndef FFT_HPP
#define FFT_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

/*!
 * \brief fft One-dimensional DFT for matrices.
 * Transform each column of the complex matrix
 * using a discrete fast Fourier transform.
 * \param X A complex matrix.
 * \return A complex matrix, whose columns are transformed
 * under DFT.
 */
Eigen::MatrixXcd fft(const Eigen::MatrixXcd& X) {
  const long m = X.rows(), n = X.cols();

  Eigen::MatrixXcd Y(m, n);

  Eigen::FFT<double> fft;

  for (long j = 0; j < n; ++j) {
    Eigen::VectorXcd Xj = X.col(j);
    Y.col(j) = fft.fwd(Xj);
  }
  return Y;
}

/*!
 * \brief fft One-dimensional inverse DFT for matrices.
 * Transform each column of the complex matrix
 * using a inverse discrete fast Fourier transform.
 * \param X A complex matrix.
 * \return A complex matrix, whose columns are transformed
 * under inverse DFT.
 */
Eigen::MatrixXcd ifft(const Eigen::MatrixXcd& X) {
  const long m = X.rows(), n = X.cols();

  Eigen::MatrixXcd Y(m, n);

  Eigen::FFT<double> fft;

  for (long j = 0; j < n; ++j) {
    Eigen::VectorXcd Xj = X.col(j);
    Y.col(j) = fft.inv(Xj);
  }
  return Y;
}

/*!
 * \brief fft Two-dimensional DFT for matrices.
 * Performs 2-dimensional DFT, using one-dimensional
 * DFT.
 * \param X A complex matrix.
 * \return A complex matrix, with Fourier coeffficients of X.
 */
Eigen::MatrixXcd fft2(const Eigen::MatrixXcd& X) {
  return fft(fft(X).transpose()).transpose();
}

/*!
 * \brief fft Two-dimensional inverse DFT for matrices.
 * Performs 2-dimensional inverse DFT,
 * using one-dimensional DFTs.
 * \param X A complex matrix.
 * \return A complex matrix, with Fourier coeffficients of X.
 */
Eigen::MatrixXcd ifft2(const Eigen::MatrixXcd& X) {
  return ifft(ifft(X).transpose()).transpose();
}

/*!
 * \brief fftr One-dimensional DFT for matrices.
 * Transform each column of the real matrix
 * using a discrete fast Fourier transform.
 * \param X A real matrix.
 * \return A complex matrix, whose columns are transformed
 * under inverse DFT.
 */
Eigen::MatrixXcd fftr(const Eigen::MatrixXd& X) {
  const long m = X.rows(), n = X.cols();

  Eigen::MatrixXcd Y(m, n);

  Eigen::FFT<double> fft;

  for (long j = 0; j < n; ++j) {
    Eigen::VectorXd Xj = X.col(j);
    Y.col(j) = fft.fwd(Xj);
  }

  return Y;
}

/*!
 * \brief ifftr One-dimensional inverse DFT for matrices.
 * Transform each column of the complex matrix
 * using a inverse discrete fast Fourier transform.
 * \param X A complex matrix.
 * \return A complex matrix, whose columns are transformed
 * under inverse DFT.
 */
Eigen::MatrixXd ifftr(const Eigen::MatrixXcd& X) {
  const long m = X.rows(), n = X.cols();

  Eigen::MatrixXd Y(m, n);

  Eigen::FFT<double> fft;

  for (long j = 0; j < n; ++j) {
    Eigen::VectorXcd Xj = X.col(j);
    Y.col(j) = fft.inv(Xj);
  }

  return Y;
}

/*!
 * \brief fft Two-dimensional real DFT for matrices.
 * Performs 2-dimensional real DFT, using one-dimensional
 * real/complex DFT.
 * \param X A real matrix.
 * \return A complex matrix, with Fourier coeffficients of X.
 */
Eigen::MatrixXcd fft2r(const Eigen::MatrixXd& X) {
  return fft(fftr(X).transpose()).transpose();
}

/*!
 * \brief fft Two-dimensional inverse real DFT for matrices.
 * Performs 2-dimensional inverse real DFT,
 * using one-dimensional real/complex DFTs.
 * \param X A complex matrix, assumed to be result of
 * real DFT (hermitian).
 * \return A real matrix, with Fourier coeffficients of X.
 */
Eigen::MatrixXd ifft2r(const Eigen::MatrixXcd& X) {
  return ifftr(ifft(X).transpose()).transpose();
}

#endif
