#include <Eigen/Dense>
#include <iostream>

namespace blur {


using Eigen::MatrixXd;

/* SAM_LISTING_BEGIN_0 */
inline MatrixXd blur(const MatrixXd &P, const MatrixXd &S) {
  typedef Eigen::Index index_t;
  auto dimensions = std::make_tuple(P.rows(), P.cols(), S.rows(), S.cols());
  const auto [ m, n, M, N ] = dimensions; 
  const index_t L = (M - 1) / 2;

  if (M != N) {
    std::cout << "Error: S not quadratic!\n";
  }

  MatrixXd C(m, n);
  for (index_t l = 1; l <= m; ++l) {
    for (index_t j = 1; j <= n; ++j) {
      double s = 0;
      for (index_t k = 1; k <= (2 * L + 1); ++k) {
        for (index_t q = 1; q <= (2 * L + 1); ++q) {
          index_t kl = l + k - L - 1;
          if (kl < 1) {
            kl += m;
          } else if (kl > m) {
            kl -= m;
          }
          index_t jm = j + q - L - 1;
          if (jm < 1) {
            jm += n;
          } else if (jm > n) {
            jm -= n;
          }
          s += P(kl - 1, jm - 1) * S(k - 1, q - 1);
        }
      }
      C(l - 1, j - 1) = s;
    }
  }
  return C;
}
/* SAM_LISTING_END_0 */


} //namespace blur
