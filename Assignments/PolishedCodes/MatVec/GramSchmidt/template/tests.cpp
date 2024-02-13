#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include "copy.hpp"
#include "doctest.h"

std::string highlight_diff(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B, std::string color_code)
{
  assert(A.rows() == B.rows() && A.cols() == B.cols());

  std::stringstream ss;
  ss.precision(4);
  ss << std::fixed;
  std::string sep = " ";
  for (int i = 0; i < A.rows(); i++)
  {
    for (int j = 0; j < A.cols(); j++)
    {
      bool is_eps_close = std::abs(A(i, j) - B(i, j)) < 1e-6;
      const std::string color = is_eps_close ? "" : color_code;
      const std::string del_color = is_eps_close ? "" : "\033[0m";
      ss << color << std::right << std::setw(7) << A(i, j) << del_color << sep;
    }
    ss << std::endl;
  }

  return ss.str();
}

void run_subcases(const Eigen::MatrixXd A)
{
  const std::string red_bgd_code = "\033[41m"; // red background color
  const std::string red_txt_code = "\033[31m"; // red text color
  const std::string reset_c_code = "\033[0m";  // reset color
  Eigen::MatrixXd Q_sol, Q_stud;
  Q_sol = gram_schmidt(A);
  Q_stud = gram_schmidt_TEST(A);

  SUBCASE("Matrix dimension checks")
  {
    std::string matrix_dim_msg = red_txt_code + "The matrix dimensions are wrong" + reset_c_code;
    REQUIRE_MESSAGE(Q_sol.rows() == Q_stud.rows(), matrix_dim_msg);
    REQUIRE_MESSAGE(Q_sol.cols() == Q_stud.cols(), matrix_dim_msg);
  }

  SUBCASE("Column orthogonality check")
  {
    for (int i = 0; i < Q_stud.cols(); ++i)
    {
      for (int j = i + 1; j < Q_stud.cols(); ++j)
      {
        std::string col_orth_msg = red_txt_code + "Columns " + std::to_string(i) + " and " + std::to_string(j) + " of Q are not orthogonal. " + reset_c_code;
        CHECK_MESSAGE(Q_stud.col(i).transpose() * Q_stud.col(j) < 1e-6, col_orth_msg);
      }
    }
  }

  SUBCASE("Column norm check")
  {
    for (int i = 0; i < Q_stud.cols(); ++i)
    {
      std::string col_norm_msg = red_txt_code + "Columns " + std::to_string(i) + " of Q is not normalized. " + reset_c_code;
      CHECK_MESSAGE(Q_stud.col(i).norm() == doctest::Approx(1.0).epsilon(1e-6), col_norm_msg);
    }
  }

  SUBCASE("Frobenius norm checks")
  {
    std::string frobenius_norm_msg = red_txt_code + "Student solution matrix different from master solution as indicated by Frobenius norm. Wrong matrix entries are marked in " + reset_c_code + red_bgd_code + "RED" + reset_c_code + ". \n";
    std::string high_lighted_diffs = highlight_diff(Q_stud, Q_sol, red_bgd_code);
    frobenius_norm_msg.append(high_lighted_diffs);
    CHECK_MESSAGE((Q_sol - Q_stud).norm() == doctest::Approx(0.).epsilon(1e-6), frobenius_norm_msg);
  }
}

TEST_SUITE("Gram-Schmidt")
{
  TEST_CASE("Eigen::MatrixXd gram_schmidt" *
            doctest::description("gram_schmidt(A)"))
  {
    SUBCASE("Matrix A")
    {
      Eigen::Matrix3d A = Eigen::Matrix3d::Identity();
      A(2, 1) = -5.0;
      A(0, 2) = 1.0;
      run_subcases(A);
    }

    SUBCASE("Matrix B")
    {
      const Eigen::MatrixXd B = Eigen::MatrixXd::Random(5, 5);
      run_subcases(B);
    }

    SUBCASE("Matrix C")
    {
      const Eigen::MatrixXd C = Eigen::MatrixXd::Random(8, 8);
      run_subcases(C);
    }
  }

  TEST_CASE("bool testGramSchmidt" * doctest::description("Test by student"))
  {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
