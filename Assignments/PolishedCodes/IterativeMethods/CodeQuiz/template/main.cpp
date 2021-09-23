///
/// Minimal runner for (9-2) Code Quiz
///
#include <iomanip>
#include <iostream>

#include "codequiz.hpp"

int main(void) {
  // You can test myfunction() for different input values:
  std::cout << std::setprecision(15) << "myfunction(10) = " << myfunction(10)
            << std::endl;

  std::cout << std::setprecision(15)
            << "myfunction_modified(10) = " << myfunction_modified(10)
            << std::endl;
}
