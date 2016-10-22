# include <iostream>
# include <Eigen/Dense>
using Eigen::MatrixXd;

/* SAM_LISTING_BEGIN_0 */
MatrixXd blur(const MatrixXd& P, const MatrixXd& S) {
  const long m = P.rows(), n = P.cols(),
             M = S.rows(), N = S.cols(), L = (M-1)/2;
  if (M != N) {
    std::cout << "Error: S not quadratic!\n";
  }

  MatrixXd C(m,n);
  for (long l = 0; l < m; ++l) {
    for (long j = 0; j < n; ++j) {
      double s = 0;
      for (long k = 0; k < (2*L+1); ++k) {
        for (long q = 0; q < (2*L+1); ++q) {
          double kl = l + k - L - 1;
          if (kl < 0) kl += m;
          else if (kl >= m) kl -= m;
          double jm = j + q - L - 1;
          if (jm < 0) jm += n;
          else if (jm >= n) jm -= n;
          s += P(kl, jm)*S(k,q);
        }
      }
      C(l,j) = s;
    }
  }
  return C;
}
/* SAM_LISTING_END_0 */
