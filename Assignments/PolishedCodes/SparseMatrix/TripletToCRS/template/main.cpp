#include <iostream>

#include "triplettoCRS.hpp"

/**
 * @brief Tests if the TripletMatrix to CRSMatrix conversion method works.
 *
 * @param n dimension of the matrix (n by n)
 * @return true if the method works
 * @return false otherwise
 */
/* SAM_LISTING_BEGIN_5 */
bool testTripletToCRS(std::size_t n) {
  bool areTheSame = false;
  // TODO: (2-13.d) Test if the conversion method from above works. To that end,
  // you may use the densify() methods.
  // START

  // END
  return areTheSame;
}
/* SAM_LISTING_END_5 */

int main() {
  constexpr std::size_t n = 1;
  const bool ret = testTripletToCRS(n);
  if (ret) {
    std::cout << "testTripletToCRS() returns true." << std::endl;
  } else {
    std::cout << "testTripletToCRS() returns false." << std::endl;
  }
}