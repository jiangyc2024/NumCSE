///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2022 SAM, D-MATH
/// Author(s): R. Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <complex>
#include <iostream>
#include <type_traits>
#include <utility>
#include <vector>

/* Related to Review Question qi:qf2 in Chapter 7 */
template <typename QUADRULE, typename FUNCTION>
auto integrate(const QUADRULE &qr, FUNCTION &&f)
    -> decltype(qr[0].second * f(qr[0].first)) {
  using Scalar = decltype(qr[0].second * f(qr[0].first));
  Scalar s{0};
  for (auto nw : qr) {
    s += nw.second * f(nw.first);
  }
  return s;
}

int main(int /*argc*/, char ** /*argv*/) {
  const double sqrt3 = std::sqrt(3) / 3.0;
  std::vector<std::pair<double, double>> twoPtsGaussRule{{-sqrt3, 1.0},
                                                         {sqrt3, 1.0}};
  std::cout << "I = "
            << integrate(twoPtsGaussRule,
                         [](double t) -> std::complex<double> {
                           return std::complex<double>(0, 1) * t * t;
                         })
            << std::endl;
  return 0;
}

/*template <typename QuadRule, typename Function>
auto  integrate(const QuadRule &qr, Function &&f) -> typename
std::result_of<Function(typename QuadRule::value_type::first_type)>::type
*/
