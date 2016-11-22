///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
// Auxiliary functions for dealing with polynomials, see \cref{sec:pi}
# include "feval.hpp"
# include "polyfit.hpp"
# include "polyval.hpp"

using Eigen::VectorXd;

int main() {
/* SAM_LISTING_BEGIN_0 */
  // Note: ``quick \& dirty'' implementation \textbf{!}
  // Lambda function representing \Blue{$x\mapsto (1+x^2)^-1$}
  auto f = [](double x) { return 1./(1 + x*x); };
  // sampling points for approximate maximum norm
  const VectorXd x = VectorXd::LinSpaced(1000, -5, 5);
  // Sample function
  const VectorXd fx = feval(f, x); // evaluate f at x

  std::vector<double> err; // Accumulate error norms here
  for (int d = 1; d <= 20; ++d) {
    // Interpolation nodes
    const VectorXd t = Eigen::VectorXd::LinSpaced(d + 1, -5, 5);
    // Interpolation data values
    const VectorXd ft = feval(f, t);
    // Compute interpolating polynomial
    const VectorXd p = polyfit(t, ft, d);
    // Evaluate polynomial interpolant
    const VectorXd y = polyval(p, x);
    // Approximate supremum norm of interpolation error
    err.push_back( (y - fx).cwiseAbs().maxCoeff() );
  }
/* SAM_LISTING_END_0 */

  mgl::Figure fig;
  fig.setlog(false, true);
  fig.plot(err, "r--+");

  fig.xlabel("degree d");
  fig.ylabel("interpolation error (maximum norm)");
  fig.save("rungeerrmax.eps");

  return 0;
}
