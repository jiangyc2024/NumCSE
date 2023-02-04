/**
 * \file main.cpp
 * \brief NPDE homework StabRK3 code
 * \author Oliver Rietmann, Philippe Peter
 * \date 13.04.2021
 * \copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <iomanip>
#include <iostream>
#include <vector>

#include "stabrk3.hpp"

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

int main() {
  // aproximate reference solution used in the convergence study
  constexpr double T = 1.0;
  Eigen::Vector2d y0(100.0, T);
  Eigen::Vector2d yT_reference = PredPrey(y0, T, std::pow(2, 14));
  std::cout << "Solution Computed by PredPrey(): "
            << yT_reference.transpose().format(CSVFormat) << std::endl;

  // Convergence study
  SimulatePredPrey();

  return 0;
}
