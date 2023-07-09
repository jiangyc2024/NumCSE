///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <cmath>

/* SAM_LISTING_BEGIN_0 */
// C++ template for simplified Newton method
template<typename Func, typename Jac, typename Vec>
void simpnewton(Vec& x, Func F, Jac DF, double rtol, double atol)
{
  auto lu = DF(x).lu(); // do LU decomposition \com{once}!
  Vec s;                // Newton correction
  double ns = NAN;      // auxiliary variables for termination control
  double nx = NAN;         
  do {
    s = lu.solve(F(x));
    x = x-s;            // new iterate
    ns = s.norm(); nx = x.norm();
  }
  // \com{termination} based on relative and absolute tolerance
  while((ns > rtol*nx) && (ns > atol));
}
/* SAM_LISTING_END_0 */
