
#include "matplotlibcpp.h"
#include "trigipcardfn.hpp"

namespace plt = matplotlibcpp;

/*!
 * \brief plot_basis Plot the shifted basis polynomials.
 * \param n $2*n+1$ will be the number of basis polynomials.
 */
/* SAM_LISTING_BEGIN_0 */

void plot_basis(int n) {
  // mesh size
  const int M = 1e3;

  // START

  // END
}
/* SAM_LISTING_END_0 */

/*!
 * \brief plot_lam Plot the Lebesgue constant $\lambda(n)$ in function of n.
 * \param points are $n = 2^k$, for k = 2,3,...,8
 * \param lambda are the Lebesgue constants.
 */

void plot_lam(std::vector<unsigned int> &points, std::vector<float> &lambda) {
  // plot using matplotlibcpp
  plt::figure();
  plt::title(" Lebesgue constant for trigonomatric interpolation ");
  plt::xlabel("n");
  plt::ylabel("λ(n)");
  plt::plot(points, lambda);
  plt::savefig("cx_out/lebesgue.png");
}
