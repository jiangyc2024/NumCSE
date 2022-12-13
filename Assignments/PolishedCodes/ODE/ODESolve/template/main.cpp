////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////

#include "odesolve.hpp"

int main() {
  auto Psi = [](double h, const Eigen::VectorXd &y0) -> Eigen::VectorXd {
    return (1 + h) * y0;
  };

  Eigen::VectorXd y0(1);
  y0 << 1.0;
  constexpr double T = 1.0;
  // Test psitilde
  Eigen::VectorXd y1 = psitilde(Psi, 1, 0.1, y0);
  std::cout << "Test psitilde\n" << y1 << std::endl;

  std::cout << "Enter 1 to test equidistant integration\n"
            << "Enter 2 to test adaptive integration\n"
            << "Enter 0 to exit\n";
  unsigned int flag;
  std::cin >> flag;
  switch (flag) {
    case 1: {
      constexpr unsigned int N = 8;

      std::cout << "Test equidistant integration" << std::endl;
      std::vector<Eigen::VectorXd> Y1 = odeintequi(Psi, T, y0, N);
      for (unsigned int i = 0; i < N; ++i) {
        std::cout << Y1[i](0) << std::endl;
      }
      double rate = testcvpExtrapolatedEuler();
      std::cout << "\nRate = " << rate << std::endl;
      break;
    }

    case 2: {
      std::cout << "\n\nTest adaptive integration" << std::endl;
      std::vector<Eigen::VectorXd> mysol =
          odeintssctrl(Psi, T, y0, 0.01, 1, 10e-5, 10e-5, 10e-5).second;
      for (unsigned int i = 0; i < 8; ++i) {
        std::cout << mysol[i](0) << std::endl;
      }

      solveTangentIVP();
      break;
    }
    default:
      break;
  }

  return 0;
}
