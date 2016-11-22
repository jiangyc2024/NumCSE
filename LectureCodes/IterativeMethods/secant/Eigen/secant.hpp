///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
/* SAM_LISTING_BEGIN_0 */
// Secand method for solving \Blue{$F(x)=0$} for \Blue{$F:D\subset\bbR\to\bbR$},
// initial guesses \Blue{$x_0,x_1$},
// tolerances \texttt{atol} (absolute), \texttt{rtol} (relative)
template<typename Func>
double secant(double x0, double x1, Func&& F, double rtol, double atol, unsigned int maxIt)
{
  double fo = F(x0);
  for (unsigned int i=0; i<maxIt; ++i) {
      double fn = F(x1);
      double s = fn*(x1-x0)/(fn-fo); // secant correction
      x0 = x1; x1 = x1-s;
      // correction based termination (relative and absolute)
      if (std::abs(s) < std::max(atol,rtol*std::min(std::abs(x0),std::abs(x1))))
	return x1;
      fo = fn; 
  }
  return x1;
}
/* SAM_LISTING_END_0 */
