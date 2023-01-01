#include <iostream>

#include "gsit.hpp"

int main() {
  constexpr unsigned int n = 9;
  const double res_norm = testGSIt(n);
  std::cout << "--> Test GSIt" << std::endl;
  std::cout << "Residual norm = " << res_norm << std::endl;
}