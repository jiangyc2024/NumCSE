#include <Eigen/Dense>
#include <iostream>
using Eigen::MatrixXd;

/* SAM_LISTING_BEGIN_0 */
MatrixXd blur(const MatrixXd &P, const MatrixXd &S) {
  const int m = P.rows(), n = P.cols(), M = S.rows(), N = S.cols(),
             L = (M - 1) / 2;
  if (M != N) {
    std::cout << "Error: S not quadratic!\n";
  }

  MatrixXd C(m, n);
  for (int l = 1; l <= m; ++l) {
    for (int j = 1; j <= n; ++j) {
      double s = 0;
      for (int k = 1; k <= (2 * L + 1); ++k) {
        for (int q = 1; q <= (2 * L + 1); ++q) {
          int kl = l + k - L - 1;
          if (kl < 1)
            kl += m;
          else if (kl > m)
            kl -= m;
          int jm = j + q - L - 1;
          if (jm < 1)
            jm += n;
          else if (jm > n)
            jm -= n;
          s += P(kl - 1, jm - 1) * S(k - 1, q - 1);
        }
      }
      C(l - 1, j - 1) = s;
    }
  }
  return C;
}
/* SAM_LISTING_END_0 */
