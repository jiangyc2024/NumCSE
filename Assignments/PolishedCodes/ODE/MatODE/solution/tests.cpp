#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>
#include <vector>

using namespace Eigen;

struct TestData {
  TestData() {
    n = 3;

    h = {0.05, 0.1, 0.3, 0.5, 1.0};

    A.resize(n, n);
    A << 0, 1, 0, 1, 0, 1, 1, 1, 0;

    Y0 = Eigen::MatrixXd::Identity(n, n);
  }

  int n;
  std::vector<double> h;
  Eigen::MatrixXd A;
  Eigen::MatrixXd Y0;
};

TestData data;

TEST_SUITE("MatODE") {
  TEST_CASE("Eigen::MatrixXd eeulstep" * doctest::description("explicit eul")) {
    for (int i = 0; i < data.h.size(); i++) {
      Eigen::MatrixXd sol = eeulstep(data.A, data.Y0, data.h[i]);
      Eigen::MatrixXd stud = eeulstep_TEST(data.A, data.Y0, data.h[i]);

      CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
    }
  }

  TEST_CASE("Eigen::MatrixXd ieulstep" * doctest::description("implicit eul")) {
    for (int i = 0; i < data.h.size(); i++) {
      Eigen::MatrixXd sol = ieulstep(data.A, data.Y0, data.h[i]);
      Eigen::MatrixXd stud = ieulstep_TEST(data.A, data.Y0, data.h[i]);

      CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
    }
  }

  TEST_CASE("Eigen::MatrixXd impstep" *
            doctest::description("implicit midpoint")) {
    for (int i = 0; i < data.h.size(); i++) {
      Eigen::MatrixXd sol = impstep(data.A, data.Y0, data.h[i]);
      Eigen::MatrixXd stud = impstep_TEST(data.A, data.Y0, data.h[i]);

      CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
    }
  }

  TEST_CASE("std::tuple<double,double,double> checkOrthogonality" *
            doctest::description("orthogonal")) {
    std::tuple<double, double, double> sol = checkOrthogonality();
    std::tuple<double, double, double> stud = checkOrthogonality_TEST();

    CHECK(std::abs(std::get<0>(sol) - std::get<0>(stud)) +
              std::abs(std::get<1>(sol) - std::get<1>(stud)) +
              std::abs(std::get<2>(sol) - std::get<2>(stud)) ==
          doctest::Approx(0.).epsilon(1e-6));
  }
}
