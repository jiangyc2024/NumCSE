#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

struct TestData {
  TestData() {
    n = 1000;
    alpha = -4.;
    beta = 1.;
    T.resize(n, n);
    std::vector<Eigen::Triplet<double> > triplets;
    triplets.reserve(3 * n - 2);
    for (std::size_t i = 0; i < n - 1; ++i) {
      triplets.push_back(Eigen::Triplet<double>(i, i, alpha));
      triplets.push_back(Eigen::Triplet<double>(i + 1, i, beta));
      triplets.push_back(Eigen::Triplet<double>(i, i + 1, beta));
    }
    triplets.push_back(Eigen::Triplet<double>(n - 1, n - 1, alpha));
    T.setFromTriplets(triplets.begin(), triplets.end());
    T.makeCompressed();
    z = Eigen::VectorXd::LinSpaced(n, -1., 1.);
    std::srand(39);
    c = T * z + 0.5 * Eigen::VectorXd::Random(n);
  }

  std::size_t n;
  double alpha, beta;
  Eigen::SparseMatrix<double> T;
  Eigen::VectorXd z, c;
};

TestData data;

TEST_SUITE("TridiagLeastSquares") {
  TEST_CASE("VectorXd lsqEst" * doctest::description("lsqEst()")) {
    Eigen::VectorXd sol = lsqEst(data.z, data.c);
    Eigen::VectorXd stud = lsqEst_TEST(data.z, data.c);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }
}
