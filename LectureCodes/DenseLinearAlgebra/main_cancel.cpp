// **********************************************************************
// Codes for NumCSE: computational costs of numerical linear algebra
// operations.
// **********************************************************************

#include <assert.h>
#include <cmath>
#include <iostream>

/* SAM_LISTING_BEGIN_0 */
double f1(double x) {
  return std::log(std::sqrt(x * x + 1) - x);
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_A */
double f1c(double x) {
  return std::log(1.0 / (std::sqrt(x * x + 1) + x));
}
/* SAM_LISTING_END_A */

/* SAM_LISTING_BEGIN_1 */
double f2(double x) {
  assert(x > 0);
  return std::log(x * x + 1) - 2 * std::log(x);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_B */
double f2c(double x) {
  assert(x > 0);
  const double y = 1.0 / x;
  return std::log(y * y + 1);
}
/* SAM_LISTING_END_B */

/* SAM_LISTING_BEGIN_2 */
double f3(double x) {
  assert((x >= -1) && (x <= 1));
  return 1 - std::sqrt(1 - x * x);
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_C */
double f3c(double x) {
  assert((x >= -1) && (x <= 1));
  const double s = x * x;
  return s / (1 + std::sqrt(1 - s));
}
/* SAM_LISTING_END_C */

/* SAM_LISTING_BEGIN_3 */
double f4(double x) {
  const double s = std::cos(x);
  return std::sqrt(1 - s * s);
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_D */
double f4c(double x) {
  return std::abs(std::sin(x));
}
/* SAM_LISTING_END_D */

int main(int argc, char **argv) {
  std::cout << "f1(1.0)= " << f1(1.0) << "<-> f1c(1.0)= " << f1c(1.0)
            << std::endl;
  std::cout << "f2(1.0)= " << f2(1.0) << "<-> f2c(1.0)= " << f2c(1.0)
            << std::endl;
  std::cout << "f3(1.0)= " << f3(1.0) << "<-> f3c(1.0)= " << f3c(1.0)
            << std::endl;
  std::cout << "f4(1.0)= " << f4(1.0) << "<-> f4c(1.0)= " << f4c(1.0)
            << std::endl;
  return 0;
}
