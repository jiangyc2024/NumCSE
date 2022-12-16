#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <vector>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    T = 1.;
    N = 1;
    y0 << 0.1, 0.2, 0.4;
    f = [](const Eigen::Vector3d &y) -> Eigen::Vector3d {
      Eigen::Vector3d fy;
      fy << y(0) * y(1), y(1) * y(2), y(2) - y(0);
      return fy;
    };
    Jf = [](const Eigen::Vector3d &y) {
      Eigen::Matrix3d J;
      J << y(1), y(0), 0, 0, y(2), y(1), -1, 0, 1;
      return J;
    };
  }

  double T;
  unsigned int N;
  Eigen::Vector3d y0;
  std::function<Eigen::Vector3d(const Eigen::Vector3d &)> f;
  std::function<Eigen::Matrix3d(const Eigen::Vector3d &)> Jf;
};

TestData data;

TEST_SUITE("CrossProduct") {
  TEST_CASE("std::vector<Eigen::VectorXd> solve_lin_mid" *
            doctest::description("Implicit linear midpoint")) {
    std::vector<Eigen::VectorXd> impl_lin_mid_sol =
        solve_lin_mid(data.f, data.Jf, data.T, data.y0, data.N);
    std::vector<Eigen::VectorXd> impl_lin_mid_stud =
        solve_lin_mid_TEST(data.f, data.Jf, data.T, data.y0, data.N);

    REQUIRE(!impl_lin_mid_stud.empty());
    REQUIRE(!impl_lin_mid_sol.empty());

    auto sol = impl_lin_mid_sol.back();
    auto stud = impl_lin_mid_stud.back();

    const bool samesize = sol.size() == stud.size();
    REQUIRE(samesize);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("void tab_crossprod" * doctest::description("Tabulate results")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
