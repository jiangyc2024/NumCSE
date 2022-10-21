#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("Sparse Approximate Inverse (SPAI)") {
  TEST_CASE("Eigen::SparseMatrix<double> spai" *
            doctest::description("Sparse Approximate Inverse")) {
    // test with two different sized random matrices
    for (std::size_t n = 20; n <= 30; n += 10) {
      Eigen::SparseMatrix<double> I(n, n);
      I.setIdentity();
      Eigen::SparseMatrix<double> M(n * n, n * n);
      Eigen::MatrixXd R = Eigen::MatrixXd::Random(n, n);
      M = Eigen::kroneckerProduct(R, I);

      Eigen::SparseMatrix<double> sol = spai(M);
      Eigen::SparseMatrix<double> stud = spai_TEST(M);

      REQUIRE(sol.rows() == stud.rows());
      REQUIRE(sol.cols() == stud.cols());
      CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-9));
    }
  }

  TEST_CASE("Eigen::SparseMatrix<double> init_A" *
            doctest::description("skipped") * doctest::skip()) {}

  // clang-format off
  TEST_CASE("std::vector<std::pair<unsigned int, unsigned int>> testSPAIPrecCG" *
      doctest::description("Iteration test")) {
    // clang-format on
    std::vector<std::pair<unsigned int, unsigned int>> sol = testSPAIPrecCG(5);
    std::vector<std::pair<unsigned int, unsigned int>> stud =
        testSPAIPrecCG_TEST(5);

    REQUIRE(sol.size() == stud.size());
    for (std::size_t i = 0; i < sol.size(); ++i) {
      CHECK(std::get<0>(sol[i]) == std::get<0>(stud[i]));
      CHECK(std::get<1>(sol[i]) == std::get<1>(stud[i]));
    }
  }
}
