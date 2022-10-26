#include <cmath>
#include <iostream>
#include <vector>

#include "chebpolyproperties.hpp"

int main() {
  constexpr unsigned int n = 10;
  std::cout << "Function checkDiscreteOrthogonality returns for n = " << n
            << ":" << std::endl
            << checkDiscreteOrthogonality(n) << std::endl;

  std::cout << "Running testBestPolyChebNodes..." << std::endl;
  testBestPolyChebNodes(n);
}