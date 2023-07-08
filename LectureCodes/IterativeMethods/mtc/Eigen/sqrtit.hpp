///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

inline
/* SAM_LISTING_BEGIN_0 */
double sqrtit(double a)
{
  double x_old = -1;
  double x = a;
  while (x_old != x) {
    x_old = x;
    x = 0.5*(x+a/x);
  }
  return x;
}
/* SAM_LISTING_END_0 */
