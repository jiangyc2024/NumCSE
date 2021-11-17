#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include "trigipcardfn.hpp"

int main() {
  // subtask 5-5.c)
  constexpr unsigned int k = 5;
  plot_basis(k);

  // subtask 5-5.h)
  constexpr unsigned int max_pt =
      6;  // the max runtime for CE,  15 originally for plots as in .pdf

  std::cout << std::setw(10) << "2^k" << std::setw(10) << "lambda(k)"
            << std::endl;

  // memory allocation for plot
  std::vector<float> lambda;
  std::vector<int> points;
  Eigen::VectorXd lam(max_pt);
  Eigen::VectorXd n(max_pt);

  for (unsigned int i = 1 << 2; i < (1 << max_pt); i <<= 1) {
    const double l = trigIpL(i);
    lambda.push_back(l);
    points.push_back(i);
    std::cout << std::setprecision(3) << std::setw(10) << i << std::setw(10)
              << l << std::endl;
  }

  // graph of the Lebesgue constant $\lambda(n)$
  plot_lam(points, lambda);
}
