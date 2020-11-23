#include "./evaldivdiff.hpp"
#include <iostream>

int main() {
  /*
   * testing divdiff code for (0) runge function and (1) sine
   */
  for (unsigned n = 4; n <= 10; n += 2) {
    // (0) : Runge function
    VectorXd t0 = VectorXd::LinSpaced(n, -5, 5),
             y0 = (1. / (1 + t0.array() * t0.array())).matrix();
    // (1) : sine
    VectorXd t1 = VectorXd::LinSpaced(n, 0, 2 * M_PI),
             y1 = t1.array().sin().matrix();
    // evaluating
    VectorXd x0 = VectorXd::LinSpaced(100, -5, 5),
             x1 = VectorXd::LinSpaced(100, 0, 2 * M_PI);
    VectorXd p0, p1;
    evaldivdiff(t0, y0, x0, p0);
    evaldivdiff(t1, y1, x1, p1);
    // output results
    std::cout << "p0 =\n " << p0.transpose() << std::endl;
    std::cout << "p1 =\n " << p1.transpose() << std::endl;
  }
  return 0;
}
