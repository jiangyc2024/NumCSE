#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Core>
#include <Eigen/SparseCore>

TEST_SUITE("SparseCCS") {
  TEST_CASE("void CCS" * doctest::description("tests val, row_ind, col_ptr")) {
    Eigen::MatrixXd A(6, 6);
    A << 4, -1, 0, -1, 0, 0, -1, 4, -1, 0, -1, 0, 0, -1, 4, 0, 0, -1, -1, 0, 0,
        4, -1, 0, 0, -1, 0, -1, 4, -1, 0, 0, -1, 0, -1, 4;
    Eigen::VectorXd val_stud, row_ind_stud, col_ptr_stud;

    double* val_sol_;
    int* row_ind_sol_;
    int* col_ptr_sol_;
    Eigen::VectorXd val_sol, row_ind_sol, col_ptr_sol;

    Eigen::SparseMatrix<double> As = A.sparseView();
    As.makeCompressed();

    val_sol_ = As.valuePtr();
    row_ind_sol_ = As.innerIndexPtr();
    col_ptr_sol_ = As.outerIndexPtr();
    val_sol = Eigen::Map<Eigen::VectorXd>(val_sol_, As.nonZeros());
    row_ind_sol =
        Eigen::Map<Eigen::VectorXi>(row_ind_sol_, As.nonZeros()).cast<double>();
    col_ptr_sol =
        Eigen::Map<Eigen::VectorXi>(col_ptr_sol_, As.cols()).cast<double>();

    CCS_TEST(A, val_stud, row_ind_stud, col_ptr_stud);

    REQUIRE(val_sol.size() == val_stud.size());
    REQUIRE(row_ind_sol.size() == row_ind_stud.size());
    REQUIRE(col_ptr_sol.size() == col_ptr_stud.size());
    CHECK((val_sol - val_stud).norm() == doctest::Approx(0.).epsilon(1e-8));
    CHECK((row_ind_sol - row_ind_stud).norm() ==
          doctest::Approx(0.).epsilon(1e-8));
    CHECK((col_ptr_sol - col_ptr_stud).norm() ==
          doctest::Approx(0.).epsilon(1e-8));
  }
}
