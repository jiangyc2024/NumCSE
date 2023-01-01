#ifndef SINH_HPP
#define SINH_HPP

#include <cmath>
#include <iomanip>
#include <iostream>

#include "sinh_stable.hpp"

/* SAM_LISTING_BEGIN_1 */
double sinh_unstable(double x) {
  const double t = std::exp(x);
  return .5 * (t - 1. / t);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void sinhError() {
  // in this solution, we also compare against an implementation based on the
  // Taylor expansion of sinh in sinh\_stable.hpp
  std::cout << "--> Relative error of various implementations:" << std::endl;
  std::cout << std::setw(10) << "k" << std::setw(15) << "own" << std::setw(15)
            << "std" << std::setw(15) << "taylor" << std::setw(15) << "err"
            << std::setw(15) << "err (taylor)" << std::endl;

  // TODO: (1-10.a) Tabulate the relative error norms for the unstable sinh
  // function. Compare against the implementation of the standard library.
  // START
  for (int k = 1; k <= 10; ++k) {
    // Value x
    const double x = std::pow<double>(10., -k);

    // Our own function
    const double mySinh = sinh_unstable(x);
    // The "standard" sinh
    const double stdSinh = std::sinh(x);
    // The stable sinh (look at the file "sinh\_stable.hpp" if
    // interested in advanced C++)
    const double taylorSinh = taylor_sinh<3>(x);

    // Relative error
    const double err = std::abs(mySinh - stdSinh) / std::abs(stdSinh);
    const double err_t = std::abs(taylorSinh - stdSinh) / std::abs(stdSinh);

    // Print error
    std::cout << std::setw(10) << k << std::setw(15) << mySinh << std::setw(15)
              << stdSinh << std::setw(15) << taylorSinh << std::setw(15) << err
              << std::setw(15) << err_t << std::endl;
  }

  // We compute and print the error bound computed using the
  // remainder of the Taylor expansion
  constexpr double x = 1e-3;
  std::cout << "--> Relative error bound for x = " << x << ":" << std::endl;
  std::cout << std::setw(10) << "n" << std::setw(15) << "error bound"
            << std::endl;
  std::cout << std::setw(10) << "1" << std::setw(15) << error_bound<1>(x)
            << std::endl;
  std::cout << std::setw(10) << "2" << std::setw(15) << error_bound<2>(x)
            << std::endl;
  std::cout << std::setw(10) << "3" << std::setw(15) << error_bound<3>(x)
            << std::endl;
  std::cout << std::setw(10) << "4" << std::setw(15) << error_bound<4>(x)
            << std::endl;
  // END
}
/* SAM_LISTING_END_2 */

#endif