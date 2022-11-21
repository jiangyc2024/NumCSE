#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>
#include <vector>

using namespace Eigen;

struct TestData {
  TestData() {
    T = 2.0;
    N = 10;

    A.resize(2, 2);
    A << 0, 0, 1, 0;

    b.resize(2);
    b << 0.5, 0.5;

    y0.resize(2);
    y0 << -1, 1;

    f = [](VectorXd y) {
      VectorXd fy(2);
      fy << -0.5 * y(0), y(0) * y(1);
      return fy;
    };
  }

  double T;
  unsigned int N;

  MatrixXd A;
  VectorXd b;
  VectorXd y0;

  std::function<VectorXd(VectorXd)> f;
};

TestData data;

TEST_SUITE("RK3Prey") {
  TEST_CASE("std::vector<State> solve" *
            doctest::description("RKIntegrator.solve()")) {
    RKIntegrator<VectorXd> sol_RK(data.A, data.b);
    RKIntegrator_TEST<VectorXd> stud_RK(data.A, data.b);

    std::vector<VectorXd> sol_vec =
        sol_RK.solve(data.f, data.T, data.y0, data.N);
    std::vector<VectorXd> stud_vec =
        stud_RK.solve_TEST(data.f, data.T, data.y0, data.N);

    bool stud_vec_exists = stud_vec.size() > 0;

    CHECK(stud_vec_exists);
    if (stud_vec_exists) {
      VectorXd sol = sol_vec.back();
      VectorXd stud = stud_vec.back();

      bool last_vec_samesize = sol.size() == stud.size();

      CHECK(last_vec_samesize);
      if (last_vec_samesize) {
        CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
      }
    }
  }

  TEST_CASE("void step" * doctest::description("RKIntegrator.step()")) {
    MESSAGE("This function wasn't tested. Rely on the solve() testing.");
  }

  TEST_CASE("double RK3prey" * doctest::description("Convergence rate")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("RKIntegrator" * doctest::description("RKIntegrator initialiser") *
            doctest::skip()) {}
}
