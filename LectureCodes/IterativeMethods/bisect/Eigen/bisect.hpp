///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
/* SAM_LISTING_BEGIN_0 */
// Searching zero of \Blue{$F$} in \Blue{$[a,b]$} by bisection
template <typename Func>
double bisect(Func&& F, double a, double b, double tol)
{
  if (a > b) std::swap(a,b); // sort interval bounds
  double fa = F(a), fb = F(b);
  if (fa*fb > 0) throw "f(a) and f(b) have same sign";
  int v=1; if (fa > 0) v=-1;
  double x = 0.5*(b+a); // determine midpoint
  // termination, relies on machine arithmetic if tol = 0
  while (b-a > tol && ((a<x) && (x<b))) // \Label[line]{bs:2}
    {
      // \Blue{$\operatorname{sgn}(f(x)) = \operatorname{sgn}(f(b))$}, then use x as next right boundary
      if (v*F(x) > 0) b=x;
      // \Blue{$\operatorname{sgn}(f(x)) = \operatorname{sgn}(f(a))$}, then use x as next left boundary
      else a=x; 
      x = 0.5*(a+b); // determine next midpoint
    }
  return x;
}
/* SAM_LISTING_END_0 */
