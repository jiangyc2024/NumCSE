#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

#include "gradientflow.hpp"

int main() {
  // Test function f:
  auto f = [](const Eigen::Vector2d &x) {
    Eigen::Vector2d y;
    y << x(0) * x(0) + x(1) * x(1), x(0) * x(1);
    return y;
  };
  // Jacobian of f:
  auto J = [](const Eigen::Vector2d &x) {
    Eigen::Matrix2d Df;
    Df << 2 * x(0), 2 * x(1), x(1), x(0);
    return Df;
  };
  // Test vector y0:
  Eigen::Vector2d y0(1, 0.5);
  // Step size for testing:
  constexpr double h = 0.5;
  // Tolerance for testing
  double tol = 1E-6;

  /*
   *  run computeStages()
   */
  std::array<Eigen::VectorXd, 5> G = computeStages(f, J, y0, h);
  std::cout << "G = computeStages(f,J,y0,h):\n";
  for (int i = 0; i < G.size(); i++) {
    std::cout << "G.at(" << i << ") = " << G.at(i).transpose() << std::endl;
  }
  /*
   *  test computeStages()
   */
  Eigen::MatrixXd Gsol(5, 2);
  Gsol << 0.234048464188, 0.0911954850442, 1.16124077723, 0.43498079856,
      0.589539233307, 0.228323917228, 0.499164120689, 0.194780654616,
      1.04447583358, 0.415950635186;
  tol = 1E-6;
  // set to false as soon as an error is encountered:
  bool passed_computeStages = true;

  for (unsigned int i = 0; i < 5; i++) {
    // Check that all vectors have been initialized
    if (G.at(i).size() == 2 && G.size() == 5) {
      // Check deviations from the correct solution
      if ((Gsol.row(i) - G.at(i).transpose()).norm() >= tol) {
        passed_computeStages = false;
        std::cout << "Test for computeStages() failed: wrong output.\n\n";
        break;
      }
    } else {
      passed_computeStages = false;
      std::cout << "Test for computeStages() failed: wrong size.\n\n";
      break;
    }
  }
  if (passed_computeStages) std::cout << "Test for computeStages() passed!\n\n";

  /*
   *  run discEvolSDIRK()
   */
  Eigen::VectorXd y1 = discEvolSDIRK(f, J, y0, h);
  std::cout << "y1 = discEvolSDIRK(f,J,y0,h) = " << y1.transpose() << std::endl;
  /*
   *  test discEvolSDIRK()
   */
  Eigen::Vector2d y1sol(2.04447583358, 0.915950635186);
  tol = 1E-6;
  if (y1.size() == y1sol.size()) {
    if ((y1 - y1sol).norm() < tol)
      std::cout << "Test for discEvolSDIRK() passed!\n\n";
    else
      std::cout << "Test for discEvolSDIRK() failed: wrong output.\n\n";
  } else
    std::cout << "Test for discEvolSDIRK() failed: wrong size.\n\n";

  /*
   *  run solveGradientFlow()
   */
  // Parameters for testing only.
  Eigen::Vector2d d;
  d << 0.5 * std::sqrt(2), 0.5 * std::sqrt(2);
  constexpr double lambda = 8.5, T = 0.5;
  constexpr unsigned int N = 10;
  y0 << 1, 0;
  std::vector<Eigen::VectorXd> Y = solveGradientFlow(d, lambda, y0, T, N);
  std::cout << "Y = solveGradientFlow(d, lambda, y0, T, N):\n";
  for (int i = 0; i < Y.size(); i++) {
    std::cout << "Y.at(" << i << ") = " << Y.at(i).transpose() << std::endl;
  }

  /*
   *  test solveGradientFlow()
   */
  tol = 1E-6;
  Eigen::MatrixXd Ysol(N + 1, 2);
  Ysol << 1, 0, 0.661574350863, -0.265183334169, 0.500522309998,
      -0.345715652452, 0.415074659032, -0.354835307771, 0.36129317252,
      -0.337895066294, 0.321645925272, -0.312568326946, 0.289168885473,
      -0.28564991388, 0.26103327209, -0.259669864129, 0.236018396524,
      -0.235490344277, 0.213530281589, -0.21332581582, 0.193223336322,
      -0.193144178833;
  // set to false as soon as an error is encountered:
  bool passed_solveGradientFlow = true;

  for (int i = 0; i <= N; i++) {
    // Check that all vectors have been initialized
    if (Y.at(i).size() == 2 && Y.size() == N + 1) {
      // Check deviations from the correct solution
      if ((Ysol.row(i) - Y.at(i).transpose()).norm() >= tol) {
        passed_solveGradientFlow = false;
        std::cout << "Test for solveGradientFlow() failed: wrong output.\n\n";
        break;
      }
    } else {
      passed_solveGradientFlow = false;
      std::cout << "Test for solveGradientFlow() failed: wrong size.\n\n";
      break;
    }
  }
  if (passed_solveGradientFlow)
    std::cout << "Test for solveGradientFlow() passed!\n\n";

  return 0;
}
