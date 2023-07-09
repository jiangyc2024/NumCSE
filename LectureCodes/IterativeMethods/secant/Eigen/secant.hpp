///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <cmath>

namespace secant {


using std::abs;
using std::max;
using std::min;

/* SAM_LISTING_BEGIN_0 */
// Secand method for solving \Blue{$F(x)=0$} for \Blue{$F:D\subset\bbR\to\bbR$},
// initial guesses \Blue{$x_0,x_1$},
// tolerances \texttt{atol} (absolute), \texttt{rtol} (relative)
template <typename Func>
double secant(double x0, double x1, Func &&F, double rtol, double atol,
              unsigned int maxIt) {
  double fo = F(x0);
  for (unsigned int i = 0; i < maxIt; ++i) {
    const double fn = F(x1);
    const double s = fn * (x1 - x0) / (fn - fo); // secant correction
    x0 = x1;
    x1 = x1 - s;
    // correction based termination (relative and absolute)
    if (abs(s) < max(atol, rtol * min(abs(x0), abs(x1)))) {
      return x1;
    }
    fo = fn;
  }
  return x1;
}
/* SAM_LISTING_END_0 */


} //namespace secant