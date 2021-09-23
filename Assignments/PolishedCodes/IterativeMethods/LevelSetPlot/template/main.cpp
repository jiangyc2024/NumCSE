///
/// Minimal runner for (9-13)
///

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "levelset.hpp"

int main() {
  auto f = [](Eigen::Vector2d x) {
    return std::pow(x(0), 2) + 2 * std::pow(x(1), 4);
  };

  /*
   * run pointLevelSet
   */
  constexpr double c = 2;
  Eigen::Vector2d x0 = {std::sqrt(2), 0};
  Eigen::Vector2d d = {2, 1};

  Eigen::Vector2d p = pointLevelSet(f, d, c, x0, 1e-10, 1e-16);
  std::cout << "pointLevelSet output:\n";
  std::cout << std::setprecision(15) << p << std::endl;

  /*
   * run pointLevelset
   */
  constexpr unsigned int n = 8;
  const double area = areaLevelSet(f, n, c);
  std::cout << "areaLevelSet output:\n";
  std::cout << area << std::endl;
}