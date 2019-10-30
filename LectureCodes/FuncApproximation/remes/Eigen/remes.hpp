#include "utils.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
// IN : f = handle to the Function, point evaluation
//      df = handle to the Derivative of f, point evaluation
//      a, b = interval boundaries
//      d = degree of polynomial
//      c will be saved to save the coefficients of the interpolant in monomial
//      basis
template <class Function, class Derivative>
void remez(const Function &f, const Derivative &df, const double a,
           const double b, const unsigned d, const double tol, VectorXd &c) {
  const unsigned n = 8 * d;                     // number of sampling points
  VectorXd xtab = VectorXd::LinSpaced(n, a, b); // points of sampling grid
  VectorXd ftab = feval(f, xtab); // function values at sampling grid
  double fsupn =
    ftab.cwiseAbs().maxCoeff();   // approximate supremum norm of \Blue{$f$}
  VectorXd dftab = feval(df, xtab); // derivative values at sampling grid

  // The vector xe stores the current guess for the alternants
  // initial guess is Chebychev alternant \eqref{remez:chebalt}
  const double h = M_PI / (d + 1);
  VectorXd xe(d + 2);
  for (unsigned i = 0; i < d + 2; ++i) {
    xe(i) = (a + b) / 2. + (a - b) / 2. * std::cos(h * i);
  }

  VectorXd fxe = feval(f, xe); // f evaluated at alternants
  const unsigned maxit = 10;
  // Main iteration loop of Remez algorithm
  for (unsigned k = 0; k < maxit; ++k) {
    // Interpolation at \Blue{$d+2$} points xe with deviations
    // \Blue{$\pm\delta$} Algorithm uses monomial basis, which is \textit{not}
    // optimal
    MatrixXd V = vander(xe), A(d + 2, d + 2);
    // build Matrix A, \com{LSE}
    A.block(0, 0, d + 2, d + 1) = V.block(0, 1, d + 2, d + 1);
    for (unsigned r = 0; r < d + 2; ++r) {
      A(r, d + 1) = std::pow(-1, r);
    }

    c = A.lu().solve(fxe); // solve for coefficients of polynomial \Blue{$q$}
    VectorXd cd(
		d); // to compute monomial coefficients of derivative \Blue{$q'$}
    for (unsigned i = 0; i < d; ++i) {
      cd(i) = (d - i) * c(i);
    }

    // Find initial guesses for the inner extremes by sampling
    // track sign changes of the derivative of the approximation error
    VectorXd deltab = polyval(cd, xtab) - dftab,
      s = deltab.head(n - 1).cwiseProduct(deltab.tail(n - 1));
    VectorXd ind = findNegative(s), // ind = find(s < 0)
      xx0 = select(xtab, ind);    // approximate zeros of e'
    const unsigned nx = ind.size(); // number of approximate zeros

    if (nx < d) { // too few extrema; bail out
      std::cerr << "Too few extrema!\n";
      return;
    }

    // \com{Secant method} to determin zeros of derivative of approximation
    // error
    VectorXd F0 = polyval(cd, xx0) - feval(df, xx0);
    // initial guesses from shifting sampling points
    VectorXd xx1 = xx0 + (b - a) / (2 * n) * VectorXd::Ones(xx0.size()),
      F1 = polyval(cd, xx1) - feval(df, xx1);
    // Main loop of the secant method
    while (F1.cwiseAbs().minCoeff() > 1e-12) {
      VectorXd xx2 = xx1 - (F1.cwiseQuotient(F1 - F0)).cwiseProduct(xx1 - xx0);
      xx0 = xx1;
      xx1 = xx2;
      F0 = F1;
      F1 = polyval(cd, xx1) - feval(df, xx1);
    }

    // Determine new approximation for alternants; store in xe
    // If too many zeros of the derivative \Blue{$(f - p)'$}
    // have been found, select those, where the deviation is maximal
    if (nx == d) {
      xe = VectorXd(xx0.size() + 2);
      xe << a, xx0, b;
    } else if (nx == d + 1) {
      xe = VectorXd(xx0.size() + 1);
      if (xx0.minCoeff() - a > b - xx0.maxCoeff())
        xe << a, xx0;
      else
        xe << xx0, b;
    } else if (nx == d + 2) {
      xe = xx0;
    } else {
      VectorXd del = (polyval(c.head(d + 1), xx0) - feval(f, xx0)).cwiseAbs(),
	ind = sort_indices(del);
      xe = select(xx0, ind.tail(d + 2));
    }

    // Deviation in sampling points and approximate alternants
    fxe = feval(f, xe);
    VectorXd del(xe.size() + 2);
    del << polyval(c.head(d + 1), a * VectorXd::Ones(1)) -
      ftab(0) * VectorXd::Ones(1),
      polyval(c.head(d + 1), xe) - fxe,
      polyval(c.head(d + 1), b * VectorXd::Ones(1)) -
      ftab(ftab.size() - 1) * VectorXd::Ones(1);
    // Approximation of supremum norm of approximation error
    const double dev = del.cwiseAbs().maxCoeff();
    // \com{Termination} of Remez iteration
    if (dev < tol * fsupn) {
      break;
    }
  }
  VectorXd tmp = c.head(d + 1);
  c = tmp;
}

/* SAM_LISTING_END_0 */
