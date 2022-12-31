#ifndef CHEBYAPPROX_HPP
#define CHEBYAPPROX_HPP

#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "intpolyval.hpp"

using complex = std::complex<double>;

/**
 * @brief Approximate the length of the ellipse
 * $\gamma_p(\theta):=cos(\theta-i \log(\rho),\forall
 * 0\leq\theta\leq2\pi,\rho>1$ by sampling for N equidistant arguments \theta
 *
 * @param rho the parameter of the ellipse curve
 * @param N number of equidistant points for approximation
 * @return double Value of the length of the ellipse
 */
/* SAM_LISTING_BEGIN_0 */
double lengthEllipse(double rho, unsigned int N) {
  double length = 0.0;
  // TODO: (6-5.a) Approximate the length of ellipse given $\rho$ and N
  // equidistant points (in $\theta$)
  // START
  const complex i(0, 1);  // imaginary unit
  // Define $\gamma_p$
  auto f = [rho, i](double theta) { return cos(theta - i * log(rho)); };
  Eigen::VectorXd eval_points = Eigen::VectorXd::LinSpaced(N + 1, 0, 2 * M_PI);
  // Calculating distance between last and current point on ellipse for all
  // points
  complex last = f(eval_points(0));
  for (unsigned int i = 0; i < N; ++i) {
    complex current = f(eval_points(i + 1));
    length += std::abs(current - last);
    last = current;
  }

  // Alternative method using an approximation over the integral of the length
  // of $\gamma_p'$.
  // $df=|\gamma_p'(\theta)|$
  double length_alt = 0.0;
  auto df = [rho](long double theta) {
    return 0.5 * sqrt(rho * rho + 1 / (rho * rho) - 2 * cos(2 * theta));
  };
  double dtheta = 2 * M_PI / N;
  double theta_i = 0;
  for (unsigned int i = 0; i < N; ++i) {
    length_alt += df(theta_i);
    theta_i += dtheta;
  }
  length_alt *= dtheta;
  // END
  return length;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Returns an approximation of
 * $\min_\rho \frac{\vert\gamma_\rho\vert \max_{z \in \gamma_p} \vert f(z)
 * \vert}{\pi d([-1,1], \gamma)}$ and the $\rho$ for which the minimum is
 * attained
 *
 * @return std::pair<double, double> pair (min, argmin)
 */
/* SAM_LISTING_BEGIN_1 */
std::pair<double, double> bestBound() {
  const double maxRho = exp(asinh(M_PI / 3.));
  constexpr unsigned int N = 1000;
  double m = 0.0, r = 0.0;  // min and rho for output
  // TODO: (6-5.d) Compute (6.5.4) with d from (6.5.5) and
  // $\vert\gamma_\rho\vert$ from (6-5.a)
  // START
  // define f and gamma
  const complex i(0, 1);  // imaginary unit
  auto f = [](complex z) { return 1. / (1. + std::exp(-3. * z)); };
  auto gamma = [i](double rho, double theta) {
    return std::cos(theta - i * std::log(rho));
  };

  m = std::numeric_limits<double>::max();  // save min in here
  // discretize rhos and remove first (open interval!)
  Eigen::VectorXd rhos = Eigen::VectorXd::LinSpaced(N + 1, 1, maxRho).head(N);
  r = rhos(0);  // save argmin here

  // discretize the ellipse
  Eigen::VectorXd thetas = Eigen::VectorXd::LinSpaced(N, 0, 2 * M_PI);

  // loop over all rhos and get compute expression
  for (unsigned int k = 0; k < rhos.size(); ++k) {
    const double rho = rhos(k);  // current rho
    // get values on the ellipse
    Eigen::VectorXcd z(thetas.size());
    for (unsigned int j = 0; j < thetas.size(); ++j) {
      z(j) = gamma(rho, thetas(j));
    }

    // sup of f of $\gamma_\rho$
    const Eigen::VectorXcd fz = z.unaryExpr(f);
    const double M = lengthEllipse(rho, N) * fz.cwiseAbs().maxCoeff() *
                     M_1_PI;  // nominator of expression

    // Distance is the one between -1 and the minimal real part of the ellipse
    // (i.e. semi-major-axis - 1):
    const double d = (rho + 1. / rho) / 2. - 1;  // denominator of expression

    if (M / d < m) {
      r = rho;  // update rho
    }
    m = std::min(m, M / d);  // update minimum
  }
  // END
  return std::make_pair(m, r);
}
/* SAM_LISTING_END_1 */

/*!
 * \brief Tabulate the estimated upperbound and actual error using 10000 points
 * for LInf error of Chebychev interpolation with order from 4 to n_max
 * \param n_max the maximum order for interpolation in the table
 * \return the pair of the list of calculated error and the list of estimated
 * upperbound
 */
/* SAM_LISTING_BEGIN_2 */
std::pair<std::vector<double>, std::vector<double>> compareErrorAndBound(
    unsigned int n_max) {
  std::vector<double> upperBound, errInf;
  constexpr unsigned N = 10000;

  // TODO: (6-5.f) Tabulate the L-inf error and the bounds.
  // START
  auto m_rho = bestBound();
  const double m = m_rho.first;
  const double rho = m_rho.second;

  // Test for degrees 4 to $n_max$
  Eigen::VectorXd degrees = Eigen::VectorXd::LinSpaced(n_max - 3, 4, n_max);

  for (unsigned int j = 0; j < degrees.size(); ++j) {
    upperBound.emplace_back(m / (std::pow(rho, degrees(j) + 1) - 1));
  }

  // Compute actual error for any degree
  auto f = [](double x) { return 1. / (1 + std::exp(-3 * x)); };

  for (unsigned int l = 0; l < degrees.size(); ++l) {
    const int n = degrees(l);
    // Cheby nodes and value of f at the nodes
    Eigen::ArrayXd k = Eigen::ArrayXd::LinSpaced(n + 1, 0, n);
    Eigen::VectorXd t = ((2 * k + 1) / (2 * (n + 1)) * M_PI).cos().matrix(),
                    y = t.unaryExpr(f);
    // Evaluate error at a couple of points
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N, -1, 1),
                    F = x.unaryExpr(f);

    // Use barycentric formula (in intpolyval.hpp)
    Eigen::VectorXd LnF = intpolyval(t, y, x);
    errInf.push_back((F - LnF).cwiseAbs().maxCoeff());
  }
  // END
  return {errInf, upperBound};
}
/* SAM_LISTING_END_2 */

#endif