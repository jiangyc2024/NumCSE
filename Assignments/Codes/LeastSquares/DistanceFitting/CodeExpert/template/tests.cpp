#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

#include "timer.h"

struct TestData {
  TestData() {
    n = 4;
    D = MatrixXd(4, 4);
    D << 0.0, -3.0, -4.0, -2.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 2.0, 0.0,
        0.0, 0.0, 0.0;

    b = VectorXd(n * (n - 1) / 2);
    int curr = 0;
    for (unsigned int i = 0; i < n - 1; i++) {
      for (unsigned int k = i + 1; k < n; k++) {
        b(curr) = D(i, k);
        curr++;
      }
    }
  }

  unsigned int n;
  Eigen::MatrixXd D;
  Eigen::VectorXd b;
};

TestData data;

TEST_SUITE("DistanceFitting") {
  TEST_CASE("SparseMatrix<double> initA" * doctest::description("initA()")) {
    Eigen::MatrixXd sol, stud;

    sol = initA(data.n);
    stud = initA_TEST(data.n);

    REQUIRE(sol.size() == stud.size());

    // Test for each row of Adense if it exists in Atest
    unsigned int n = sol.rows();
    bool all_rows_contained = true;  // True if stud is a permutation of sol

    for (unsigned int i = 0; i < n; i++) {  // Loop over all rows of stud
      VectorXd row = stud.row(i).transpose();
      bool contained = false;  // True if row i is contained in sol
      for (unsigned int j = 0; j < n; j++) {  // Loop over all rows of sol
        VectorXd row_test = sol.row(j).transpose();
        if ((row - row_test).norm() == 0) contained = true;
      }
      if (contained == false) all_rows_contained = false;
    }

    CHECK(all_rows_contained);
  }
  TEST_CASE("VectorXd solveExtendedNormalEquations" *
            doctest::description("solveExtendedNormalEquations()")) {
    VectorXd sol, stud;

    stud = solveExtendedNormalEquations_TEST(data.D);

    MatrixXd A = initA(data.n);
    sol = (A.transpose() * A).lu().solve(A.transpose() * data.b);

    REQUIRE(stud.size() == sol.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-9));
  }
  TEST_CASE("VectorXd solveNormalEquations" *
            doctest::description("solveNormalEquations()")) {
    // Checking for correct solution
    VectorXd sol, stud;

    stud = solveNormalEquations_TEST(data.D);
    MatrixXd A = initA(data.n);
    sol = (A.transpose() * A).lu().solve(A.transpose() * data.b);

    REQUIRE(stud.size() == sol.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-9));

    MESSAGE(
        "This is a test for the runtime of your implementation. There might be "
        "a timeout occuring if your implementation is to slow.");
    // Checking for correct runtime
    int steps = 6;  // testing until dim = 2^(steps+2)
    VectorXd time(steps);

    for (int s = 0; s < steps; s++) {
      int dim = std::pow(2, s + 3);
      MatrixXd rand = MatrixXd::Zero(dim, dim);
      rand.triangularView<Upper>() = MatrixXd::Random(dim, dim);
      rand.diagonal() = VectorXd::Zero(dim);

      Timer timer = Timer();
      for (int i = 0; i < 10; i++) {
        timer.start();
        VectorXd y_test = solveNormalEquations(rand);
        timer.stop();
      }
      time(s) = timer.mean();
    }
    VectorXd time_log = time.array().log().matrix();
    VectorXd coeff = time_log.tail(steps - 1) - time_log.head(steps - 1);
    coeff = coeff / std::log(2);

    CHECK(coeff.maxCoeff() < 2.5);
  }
  TEST_CASE("bool testNormalEquations" *
            doctest::description("testNormalEquations()")) {
    CHECK(testNormalEquations_TEST(data.D));
  }
}
