#include "levelset.hpp"

int main() {

  auto f = [](Vector2d x) { return std::pow(x(0), 2) + 2 * std::pow(x(1), 4); };

  /*
   * run pointLevelSet
   */
  double c = 2;

  Vector2d x0 = {std::sqrt(2), 0};
  Vector2d x1 = {0, 1};
  Vector2d d = {2, 1};

  Vector2d p = pointLevelSet(f, d, c, x0, x1, 1e-10, 1e-16);
  std::cout << "pointLevelSet output:\n";
  std::cout << std::setprecision(15) << p << "\n";

  /*
   * test of pointLevelSet
   */
  Vector2d ptest = {1.28718850581117, 0.643594252905583};
  if ((p - ptest).norm() < 1e-14) {
    std::cout << "Test for pointLevelSet passed!\n\n";
  } else
    std::cout << "Test for pointLevelSet failed: wrong output!\n\n";

  /*
   * run areaLevelSet
   */

  unsigned int n = 8;
  double area = areaLevelSet(f, n, c);
  std::cout << "areaLevelSet output:\n";
  std::cout << area << "\n";

  /*
   * test of areaLevelSet
   */

  double areatest = 4.26647319712633;
  if (std::abs(area - areatest) < 1e-14) {
    std::cout << "Test for areaLevelSet passed!\n\n";
  } else
    std::cout << "Test for areaLevelSet failed: wrong output!\n\n";

}
