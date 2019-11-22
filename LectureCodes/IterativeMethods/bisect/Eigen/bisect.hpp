///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Wolf,
//             Filippo Leonardi <filippo.leonardi@sam.math.ethz.ch>,
//             Till Ehrengruber <tille@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <type_traits>
#include <cassert>

/* SAM_LISTING_BEGIN_0 */
// Searching zero of \Blue{$F$} in \Blue{$[a,b]$} by bisection
template <typename Func, typename Scalar>
Scalar bisect(Func&& F, Scalar a, Scalar b, Scalar tol)
{
  if (a > b) std::swap(a,b); // sort interval bounds
  if (F(a)*F(b) > 0) throw "f(a) and f(b) have same sign";
  static_assert(std::is_floating_point<Scalar>::value,
               "Scalar must be a floating point type");
  int v=F(a) < 0 ? 1 : -1;
  Scalar x = (a+b)/2; // determine midpoint
  // termination, relies on machine arithmetic if tol = 0
  while (b-a > tol) { // \Label[line]{bs:2}
    assert(a<=x && x<=b); // assert invariant
    // \Blue{$\operatorname{sgn}(f(x)) = \operatorname{sgn}(f(b))$}, then use x as next right boundary
    if (v*F(x) > 0) b=x;
    // \Blue{$\operatorname{sgn}(f(x)) = \operatorname{sgn}(f(a))$}, then use x as next left boundary
    else a=x; 
    x = (a+b)/2; // determine next midpoint
  }
  return x;
}
/* SAM_LISTING_END_0 */
