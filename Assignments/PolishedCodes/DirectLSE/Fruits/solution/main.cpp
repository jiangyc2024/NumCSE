#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "fruits.hpp"

int main() {
  Eigen::VectorXd prices = fruitPrice();
  // Print the result nicely
  std::vector<std::string> fruit_names{"Mango",  "Kiwi",        "Lychee",
                                       "Banana", "Pomegranate", "Pineapple"};
  std::cout << std::fixed << std::setprecision(2);  // stream formatting
  std::cout << "The fruits cost:" << std::endl;
  for (unsigned int i = 0; i < 6; ++i) {
    std::cout << fruit_names[i] << ": " << prices[i] << " sFr." << std::endl;
  }
}