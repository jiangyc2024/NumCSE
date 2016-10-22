# include <iostream>
# include <Eigen/Dense>
using Eigen::MatrixXd;

MatrixXd blur(const MatrixXd& P, const MatrixXd& S) {
  const long m = P.rows(), n = P.cols(),
             M = S.rows(), N = S.cols(), L = (M-1)/2;
  if (M != N) {
    std::cout << "Error: S not quadratic!\n";
  }

  MatrixXd C(m,n);
  for (long l = 1; l <= m; ++l) {
    for (long j = 1; j <= n; ++j) {
      double s = 0;
      for (long k = 1; k <= (2*L+1); ++k) {
        for (long q = 1; q <= (2*L+1); ++q) {
          double kl = l + k - L - 1;
          if (kl < 1) kl += m;
          else if (kl > m) kl -= m;
          double jm = j + q - L - 1;
          if (jm < 1) jm += n;
          else if (jm > n) jm -= n;
          s += P(kl-1, jm-1)*S(k-1,q-1);
        }
      }
      C(l-1,j-1) = s;
    }
  }
  return C;
}
