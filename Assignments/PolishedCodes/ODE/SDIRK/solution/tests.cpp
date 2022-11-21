#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>
#include <vector>

struct TestData {
  TestData() {
    z0 << 1.0, 0.0;
    gamma = {1.0, (3 + std::sqrt(3)) / 6.0};
    T = {5.0, 6.0, 7.0, 9.0, 10.0};
    N = {20, 40, 80, 160, 320};
  }

  Eigen::Vector2d z0;
  double h;
  std::vector<double> gamma;
  std::vector<double> T;
  std::vector<unsigned int> N;
};

TestData data;

TEST_SUITE("SDIRK") {
  TEST_CASE("Eigen::Vector2d sdirkStep" * doctest::description("sdirkStep")) {
    for (int i = 0; i < data.gamma.size(); i++) {
      for (int j = 0; j < data.T.size(); j++) {
        for (int k = 0; k < data.N.size(); k++) {
          // Use equidistant timesteps
          double h = data.T[j] / data.N[k];

          Eigen::Vector2d sol = sdirkStep(data.z0, h, data.gamma[i]);
          Eigen::Vector2d stud = sdirkStep_TEST(data.z0, h, data.gamma[i]);

          CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
        }
      }
    }
  }

  TEST_CASE("std::vector<Eigen::Vector2d> sdirkSolve" *
            doctest::description("sdirkSolve")) {
    for (int i = 0; i < data.gamma.size(); i++) {
      for (int j = 0; j < data.T.size(); j++) {
        for (int k = 0; k < data.N.size(); k++) {
          std::vector<Eigen::Vector2d> sol_vec =
              sdirkSolve(data.z0, data.N[k], data.T[j], data.gamma[i]);
          std::vector<Eigen::Vector2d> stud_vec =
              sdirkSolve_TEST(data.z0, data.N[k], data.T[j], data.gamma[i]);

          // Check final value
          Eigen::Vector2d sol = sol_vec.back();
          Eigen::Vector2d stud = stud_vec.back();

          CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
        }
      }
    }
  }

  TEST_CASE("double cvgSDIRK" * doctest::description("cvgSDIRK")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
