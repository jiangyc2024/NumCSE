///
/// Minimal runner for task (8-8.c)
///

#include <iostream>
#include "gauss_pts.hpp"

int main() {
  std::cout << "g(xi^21_11) = " << testCompGGaussPts() << std::endl;
}
