
#include<iostream>
#include "julia.hpp"

int main() {
  // Tests will only pass if z = (-0.888,1.333).
  Vector2d z(-0.888,1.333);
  std::cout << "Test: F_vector(z) = " << F_vector(z).transpose() << "\n\n";
  std::cout << "Test: DF_matrix(z) = \n" << DF_matrix(z) << "\n\n";
  
  // No tests.
  julia();
}
