/***********************************************************************
 *                                                                     *
 * Demo code                                                           *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: R.H.                                                        *
 * Date: September 2023                                                *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

// Header for basic IO
#include <iostream>
// Provides random acccess container class
#include <string>
#include <utility>

/* SAM_LISTING_BEGIN_1 */
int main() {
  // Some local variables
  double x = 1.0;
  int n = 2;
  std::string s("x = ");

  auto capture_by_copy = [x](double y) -> double {
    // x *= 2; ERROR: "Copied" variables are read only in lambda functions.
    // The behave like const references!
    return x + y; };
  auto capture_all_by_copy = [=](double y) -> double {
    std::cout << "\t >>>> " << s << x << " >>>> " << std::flush;
    return x + y; };
  auto capture_by_reference = [&x](double y) -> double {
    x *= 2; return x + y; };
  auto capture_by_const_reference = [&x =
                                         std::as_const(x)](double y) -> double {
    // x *= 2; ERROR: cannot write to const reference
    // Same as if you had simply captured x as "copy"
    return x + y;
  };
  auto capture_all_by_reference = [&](double y) -> double {
    std::cout << "\t >>>> " << s << x << " >>>> " << std::flush;
    x *= 2; n++; return x + y;  };

  std::cout << "I. "
            << " x = " << x << ", n = " << n << std::endl;
  std::cout << "Capture by copy = " << capture_by_copy(3.14) << std::endl;
  std::cout << "II. "
            << " x = " << x << ", n = " << n << std::endl;
  std::cout << "Capture all by copy = " << capture_all_by_copy(3.14)
            << std::endl;
  std::cout << "III. "
            << " x = " << x << ", n = " << n << std::endl;
  std::cout << "Capture by reference = " << capture_by_reference(3.14)
            << std::endl;
  std::cout << "IV. "
            << " x = " << x << ", n = " << n << std::endl;
  std::cout << "Capture by const reference = "
            << capture_by_const_reference(3.14) << std::endl;
  std::cout << "V. "
            << " x = " << x << ", n = " << n << std::endl;
  std::cout << "Capture all by reference = " << capture_all_by_reference(3.14)
            << std::endl;
  std::cout << "VI. "
            << " x = " << x << ", n = " << n << std::endl;
  return (0);
}
/* SAM_LISTING_END_1 */
