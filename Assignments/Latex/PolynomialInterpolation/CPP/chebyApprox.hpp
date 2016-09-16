# include <Eigen/Dense>
# include <iostream>
# include <limits>
# include <cmath>
# include <complex>

using Eigen::VectorXd;
using Eigen::VectorXcd;
typedef std::complex<double> complex;

// Returns an approximation of $\min_\rho \frac{\max_{z \in \gamma_p} \vert f(z) \vert}{d([-1,1], \gamma)}$
// and the $\rho$ for which the minimum is attained
// IN : maxRho   = maximal value of tho, compute min in (1, maxRho)
//      N        = number of intervals to discretize
//      f, gamma = function handles for f and gamma
std::pair<double, double> chebyApprox(const double maxRho, const unsigned N) {

  // define f and gamma
  const complex i(0,1); // imaginary unit
  auto f = [](complex z) { return 1./(1. + std::exp(-3.*z)); };
  auto gamma = [i](double rho, double theta) {
    return std::cos(theta - i*std::log(rho));
  };

  double m = std::numeric_limits<double>::max(); // save min in here
  // discretize rhos and remove first and last value (open interval!)
  VectorXd rhosClosedIntervall = VectorXd::LinSpaced(N, 1, maxRho),
           rhos = rhosClosedIntervall.segment(1, N-2);
  double r = rhos(0); // save argmin here 

  // discretize the ellipse
  VectorXd thetas = VectorXd::LinSpaced(N, 0, 2*M_PI);

  // loop over all rhos and get compute expression
  for (int k = 0; k < rhos.size(); ++k) {
    const double rho = rhos(k); // current rho
    // get values on the ellipse
    VectorXcd z(thetas.size());
    for (int j = 0; j < thetas.size(); ++j) {
      z(j) = gamma(rho, thetas(j));
    }

    // sup of f of $\gamma_\rho$
    const VectorXcd fz = z.unaryExpr(f);
    const double M = fz.cwiseAbs().maxCoeff(); // nominator of expression

    // Distance is the one between -1 and the minimal real part of the ellipse (i.e. semi-major-axis - 1):
    const double d = (rho + 1./rho)/2. - 1; // denominator of expression

    if (M/d < m) {
      r = rho; // update rho
    }
    m = std::min(m, M/d); // update minimum
  }

  return std::make_pair(m, r);
}
