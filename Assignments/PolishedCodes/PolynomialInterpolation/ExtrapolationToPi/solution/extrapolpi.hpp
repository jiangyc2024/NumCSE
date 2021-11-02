#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <string>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/* SAM_LISTING_BEGIN_0 */
double extrapolate_to_pi(const unsigned int k) {
  double pi = 0.0;

  // TODO (5-8.b): Use the Aitken-Neville scheme to approximate
  // pi by extrapolation.
  // Hint: You can use the constant M_PI and Eigen's sin()
  // when calculating $s_j$ for $j=2,...,k$.

  // START

  // Data for interpolation.
  // Note: We are using eigen arrays because they allow for coefficient-wise
  // operations while eigen matrices/vectors are used for linear algebra
  // operations.
  Eigen::ArrayXd J = Eigen::ArrayXd::LinSpaced(k - 1, 2, k);
  Eigen::ArrayXd T = 1.0 / J;  // $t_j = 1/j$
  Eigen::ArrayXd S = J * sin(M_PI * T);

  // Aitken-Neville scheme for evaluating the interpolating
  // polynomial of $(t_j,s_j)$ at x = 0.
  for (int l = 0; l < S.size(); l++) {
    for (int m = l - 1; m >= 0; m--) {
      S(m) = S(m + 1) + (S(m + 1) - S(m)) * (0.0 - T(l)) / (T(l) - T(m));
    }
  }
  // Now S(0) has contributions from all data points.
  pi = S(0);
  // END

  return pi;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void plotExtrapolationError(const unsigned int kmax) {
  plt::figure();
  // TODO (5-8.c): Plot the error made by extrapolate_to_pi(k) for
  // k = 2,3,...,10. Use a log-scale for the y-axis, and a linear x-axis.
  // Also, tabulate the results and errors of extrapolate_to_pi().
  // Hint: Use the constant M_PI to compute the error made by
  // extrapolate_to_pi().
  // Hint: matplotlibcpp (here = plt) implements Python's
  // matplotlib functions such as figure(), plot(), xlabel(), title(),...
  // You can use main.cpp of the LinearDataFit problem as a reference.

  // START
  Eigen::VectorXd K = Eigen::VectorXd::LinSpaced(kmax - 1, 2, kmax);
  Eigen::VectorXd Err(K.size());

  std::cout << std::setw(5) << " k" << std::setw(15) << "approx. pi"
            << std::setw(15) << "error" << std::endl;

  for (int l = 0; l < K.size(); l++) {
    double pi = extrapolate_to_pi(K(l));
    Err(l) = std::abs(pi - M_PI);
    std::cout << std::setw(5) << K(l) << std::setw(15) << pi << std::setw(15)
              << Err(l) << std::endl;
  }

  // Number of data points is k-1, so polynomial degree is k-2.
  Eigen::VectorXd Deg = K - 2 * Eigen::VectorXd::Ones(kmax - 1);
  plt::semilogy(Deg, Err);
  plt::title("Approximation of pi");
  plt::xlabel("Polynomial degree");
  plt::ylabel("Error");
  plt::grid();
  // END

  std::string path = "./cx_out/pi_error_" + std::to_string(kmax) + ".png";
  plt::savefig(path);
}
/* SAM_LISTING_END_1 */
