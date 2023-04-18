///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2023- SAM, D-MATH
/// Author(s): Ralf Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>

/* SAM_LISTING_BEGIN_1 */
double fun(double x) {
  bool neg = false;
  if (x < 0) {
    neg = true;
    x = -x;
  }
  unsigned int f = 0;
  while (x > 1) {
    f++;
    x /= 2.0;
  }
  double v = 1.0 + x;
  double num = x;
  double den = 1.0;
  for (int i = 2; true; ++i) {
    double s = (num *= x) / (den *= i);
    if (s == 0.0) {
        while (f--) {
          v = v * v;
        }
      if (neg) v = 1.0 / v;
      break;
    }
    v += s;
  }
  return v;
}
/* SAM_LISTING_END_1 */

int main(int /*argc*/, char** /*argv*/) {
  std::cout << "exp(2.0) = " << std::exp(2.0) << " <-> fun(2.0) = " << fun(2.0)
            << std::endl;
  std::cout << "exp(0.5) = " << std::exp(0.5) << " <-> fun(0.5) = " << fun(0.5)
            << std::endl;
  std::cout << "exp(-3) = " << std::exp(-3.0) << " <-> fun(-3) = " << fun(-3.0)
            << std::endl;
  return 0;
}
