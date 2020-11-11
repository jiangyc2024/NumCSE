
#include "matplotlibcpp.h"
#include "trigipcardfn.hpp"

namespace plt = matplotlibcpp;

/*!
 * @brief plot_basis Plot the shifted basis polynomials.
 * @param[in] n $2*n+1$ will be the number of basis polynomials.
 */
/* SAM_LISTING_BEGIN_0 */

void plot_basis(int n) {
  // mesh size
  const int M = 1e3;

  // TO DO: use the function trigpolyvalequid from trigipcardfn.hpp to plot the cardinal basis function $b_0(t)$ in function of $t$ for n = 5
  // START

  // basis vector $e_1$
  ArrayXd e = ArrayXd::Zero(2 * n + 1);
  e(0) = 1;
  VectorXd y;
  trigpolyvalequid(e, M, y);

  ArrayXd t = ArrayXd::LinSpaced(M, 0, 1);
  ArrayXd zer = ArrayXd::Zero(t.size());

  // Shift function right a bit
  ArrayXd y_shift(M);
  unsigned int h = M / (2 * n + 1);
  y_shift << y.tail(h), y.head(M - h);

  // plotting functions using matplotlib
  plt::figure();
  plt::title("b_0(t)");
  plt::xlabel("t");
  plt::ylabel("y");
  plt::plot(t, y_shift, "r", {{"label", "b_0(t)"}});
  plt::plot(t, zer, "--");
  plt::legend();
  plt::savefig("cx_out/basis_funtion.png");

  // END
}
/* SAM_LISTING_END_0 */

/*!
 * @brief plot_lam Plot the Lebesgue constant $\lambda(n)$ in function of n.
 * @param[in] points are $n = 2^k$, for k = 2,3,...,8
 * @param[in] lambda are the Lebesgue constants.
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
