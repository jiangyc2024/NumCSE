
#include<iostream>
#include "julia.hpp"

int main() {
  // Tests will only pass if z = (-0.888,1.333).
  Vector2d z(-0.888,1.333);
  std::cout << "Test: F(z) = " << F(z).transpose() << "\n\n";
  std::cout << "Test: DF(z) = \n" << DF(z) << "\n\n";
  
  // No tests.
  julia();
}