#include <iostream>

#include "triplettoCRS.hpp"

int main() {
  constexpr std::size_t n = 10;
  const bool ret = testTripletToCRS(n);
  if (ret) {
    std::cout << "testTripletToCRS() returns true." << std::endl;
  } else {
    std::cout << "testTripletToCRS() returns false." << std::endl;
  }
}