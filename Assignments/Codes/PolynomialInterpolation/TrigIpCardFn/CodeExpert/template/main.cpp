
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include "plot.hpp"
#include "trigipcardfn.hpp"

using namespace Eigen;

int main() {
  // subtask 6-5.c)
  const int k = 5;
  plot_basis(k);

  // subtask 6-5.h)
  const int max_pt =
      6;  // the max runntime for CE,  15 originally for plots as in .pdf

  std::cout << std::setw(10) << "2^k" << std::setw(10) << "lambda(k)"
            << std::endl;

  // memory allocation for plot
  std::vector<float> lambda;
  std::vector<int> points;
  VectorXd lam(max_pt);
  VectorXd n(max_pt);

  for (unsigned int i = 1 << 2; i < (1 << max_pt); i = i << 1) {
    double l = trigIpL(i);
    lambda.push_back(l);
    points.push_back(i);
    std::cout << std::setprecision(3) << std::setw(10) << i << std::setw(10)
              << l << std::endl;
  }

  // graph of the Lebesgue constant $\lambda(n)$
  plot_lam(points, lambda);
}
