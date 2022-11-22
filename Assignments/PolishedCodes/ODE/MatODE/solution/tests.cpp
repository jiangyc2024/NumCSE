#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <vector>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    n = 3;

    h = {0.05, 0.1, 0.3, 0.5, 1.0};

    A.resize(n, n);
    A << 0, 1, 0, 1, 0, 1, 1, 1, 0;

    Y0 = Eigen::MatrixXd::Identity(n, n);
  }

  unsigned int n;
  std::vector<double> h;
  Eigen::MatrixXd A;
  Eigen::MatrixXd Y0;
};

TestData data;

TEST_SUITE("MatODE") {
  TEST_CASE("Eigen::MatrixXd eeulstep" * doctest::description("explicit eul")) {
    for (unsigned int i = 0; i < data.h.size(); i++) {
      Eigen::MatrixXd sol = eeulstep(data.A, data.Y0, data.h[i]);
      Eigen::MatrixXd stud = eeulstep_TEST(data.A, data.Y0, data.h[i]);

      CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
    }
  }

  TEST_CASE("Eigen::MatrixXd ieulstep" * doctest::description("implicit eul")) {
    for (unsigned int i = 0; i < data.h.size(); i++) {
      Eigen::MatrixXd sol = ieulstep(data.A, data.Y0, data.h[i]);
      Eigen::MatrixXd stud = ieulstep_TEST(data.A, data.Y0, data.h[i]);

      CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
    }
  }

  TEST_CASE("Eigen::MatrixXd impstep" *
            doctest::description("implicit midpoint")) {
    for (unsigned int i = 0; i < data.h.size(); i++) {
      Eigen::MatrixXd sol = impstep(data.A, data.Y0, data.h[i]);
      Eigen::MatrixXd stud = impstep_TEST(data.A, data.Y0, data.h[i]);

      CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
    }
  }

  TEST_CASE("std::tuple<double, double, double> checkOrthogonality" *
            doctest::description("orthogonal")) {
    auto [sol_0, sol_1, sol_2] = checkOrthogonality();
    auto [stud_0, stud_1, stud_2] = checkOrthogonality_TEST();

    CHECK(std::abs(sol_0 - stud_0) + std::abs(sol_1 - stud_1) +
              std::abs(sol_2 - stud_2) ==
          doctest::Approx(0.).epsilon(1e-6));
  }
}
