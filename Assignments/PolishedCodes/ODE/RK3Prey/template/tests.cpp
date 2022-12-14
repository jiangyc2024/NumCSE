#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <vector>

#include "copy.hpp"
#include "doctest.h"

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

    f = [](Eigen::VectorXd y) {
      Eigen::VectorXd fy(2);
      fy << -0.5 * y(0), y(0) * y(1);
      return fy;
    };
  }

  double T;
  unsigned int N;

  Eigen::MatrixXd A;
  Eigen::VectorXd b;
  Eigen::VectorXd y0;

  std::function<Eigen::VectorXd(Eigen::VectorXd)> f;
};

TestData data;

TEST_SUITE("RK3Prey") {
  TEST_CASE("RKIntegrator [OUT_OF_CLASS]" *
            doctest::description("Constructor") * doctest::skip()) {}

  TEST_CASE("std::vector<State> solve [OUT_OF_CLASS]" *
            doctest::description("RKIntegrator.solve()")) {
    RKIntegrator<Eigen::VectorXd> sol_RK(data.A, data.b);
    RKIntegrator_TEST<Eigen::VectorXd> stud_RK(data.A, data.b);

    std::vector<Eigen::VectorXd> sol_vec =
        sol_RK.solve(data.f, data.T, data.y0, data.N);
    std::vector<Eigen::VectorXd> stud_vec =
        stud_RK.solve(data.f, data.T, data.y0, data.N);

    bool stud_vec_exists = stud_vec.size() > 0;

    REQUIRE(stud_vec_exists);
    Eigen::VectorXd sol = sol_vec.back();
    Eigen::VectorXd stud = stud_vec.back();

    bool last_vec_samesize = sol.size() == stud.size();

    REQUIRE(last_vec_samesize);
    if (last_vec_samesize) {
      CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
    }
  }

  TEST_CASE("void step [OUT_OF_CLASS]" *
            doctest::description("RKIntegrator.step()")) {
    MESSAGE("This function wasn't tested. Rely on the solve() testing.");
  }

  TEST_CASE("double RK3prey" * doctest::description("Convergence rate")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
