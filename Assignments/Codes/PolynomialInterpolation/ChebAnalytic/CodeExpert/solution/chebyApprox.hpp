#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <limits>
#include "intpolyval.hpp"

using Eigen::VectorXcd;
using Eigen::VectorXd;
typedef std::complex<double> complex;

/*!
 * \brief Approximate the length of the ellipse
 * $\gamma_p(\theta):=cos(\theta-i \log(\rho),\forall
 * 0\leq\theta\leq2\pi,\rho>1$ by sampling for N equidistant arguments \theta
 * \param rho, the parameter of the ellipse curve
 * \param N, number of equidistant points for approximation
 * \return Value of the length of the ellipse
 */
/* SAM_LISTING_BEGIN_0 */
double lengthEllipse(double rho, unsigned int N) {
  // TODO: (7-5.a) Approximate the length of ellipse given $\rho$ and N
  // equidistant points (in $\theta$)
  double length = 0.0;
  // START
  // $f=|\gamma_p'(\theta)|$
  auto f = [rho](long double theta) {
    return 0.5 * sqrt(rho * rho + 1 / (rho * rho) - 2 * cos(2 * theta));
  };
  double dtheta = 2 * M_PI / N;
  double theta_i = 0;
  for (unsigned int i = 0; i < N; i++) {
    length += f(theta_i);
    theta_i += dtheta;
  }
  length *= dtheta;
  // END
  return length;
}
/* SAM_LISTING_END_0 */

/*!
 * \brief Returns an approximation of
 * $\min_\rho \frac{\vert\gamma_\rho\vert \max_{z \in \gamma_p} \vert f(z)
 * \vert}{\pi d([-1,1], \gamma)}$ and the $\rho$ for which the minimum is
 * attained
 */
/* SAM_LISTING_BEGIN_1 */
std::pair<double, double> bestBound() {
  // TODO: (7-5.d) Compute (7.5.2) with d from (7.5.3) and
  // $\vert\gamma_\rho\vert$ from (7-5.a) START
  double maxRho = exp(asinh(M_PI / 3.));
  unsigned N = 1000;
  // define f and gamma
  const complex i(0, 1);  // imaginary unit
  auto f = [](complex z) { return 1. / (1. + std::exp(-3. * z)); };
  auto gamma = [i](double rho, double theta) {
    return std::cos(theta - i * std::log(rho));
  };

  double m = std::numeric_limits<double>::max();  // save min in here
  // discretize rhos and remove first (open interval!)
  VectorXd rhos = VectorXd::LinSpaced(N+1, 1, maxRho).head(N);
  double r = rhos(0);  // save argmin here

  // discretize the ellipse
  VectorXd thetas = VectorXd::LinSpaced(N, 0, 2 * M_PI);

  // loop over all rhos and get compute expression
  for (unsigned int k = 0; k < rhos.size(); ++k) {
    const double rho = rhos(k);  // current rho
    // get values on the ellipse
    VectorXcd z(thetas.size());
    for (int j = 0; j < thetas.size(); ++j) {
      z(j) = gamma(rho, thetas(j));
    }

    // sup of f of $\gamma_\rho$
    const VectorXcd fz = z.unaryExpr(f);
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
  return std::make_pair(m, r);
  // END
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
  // TODO: (7-5.f)
  // BEGIN
  const unsigned N = 10000;
  auto m_rho = bestBound();
  double m = m_rho.first;
  double rho = m_rho.second;

  // Test for degrees 4 to $n_max$
  VectorXd degrees = VectorXd::LinSpaced(n_max - 3, 4, n_max);
  std::vector<double> upperBound;

  for (int j = 0; j < degrees.size(); ++j) {
    upperBound.emplace_back(m / (std::pow(rho, degrees(j) + 1) - 1));
  }

  // Compute actual error for any degree
  auto f = [](double x) { return 1. / (1 + std::exp(-3 * x)); };
  std::vector<double> errInf;

  for (int l = 0; l < degrees.size(); ++l) {
    const int n = degrees(l);
    // Cheby nodes and value of f at the nodes
    Eigen::ArrayXd k = Eigen::ArrayXd::LinSpaced(n + 1, 0, n);
    VectorXd t = ((2 * k + 1) / (2 * (n + 1)) * M_PI).cos().matrix(),
             y = t.unaryExpr(f);
    // Evaluate error at a couple of points
    VectorXd x = VectorXd::LinSpaced(N, -1, 1), F = x.unaryExpr(f);

    // Use barycentric formula (in intpolyval.hpp)
    VectorXd LnF;
    intpolyval(t, y, x, LnF);
    errInf.push_back((F - LnF).cwiseAbs().maxCoeff());
  }
  // END
  return {errInf, upperBound};
}
/* SAM_LISTING_END_2 */
