#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>
#include <functional>

struct TestData {
  TestData() {
    f = [](double x){ return std::atan(x) - 0.123; };
    df = [](double x){ return 1. / (x * x + 1.); };
    A = Eigen::MatrixXd::Random(5, 5);
    A = (A * A.transpose()).eval();
    c = Eigen::VectorXd::Random(5);
    c = c.cwiseAbs();
  }
  
  Eigen::MatrixXd A;
  Eigen::VectorXd c;
  std::function<double(double)> f, df;
  const double a = 0.123;
};

TestData data;

TEST_SUITE("Modified Newton method") {
  TEST_CASE("double norm" * doctest::description("some helpers") * doctest::skip()) {}
  
  TEST_CASE("Scalar mod_newt_step_scalar" * doctest::description("scalar modified Newton step")) {
    const double sol = mod_newt_step_scalar(1., data.f, data.df);
    const double stud = mod_newt_step_scalar_TEST(1., data.f, data.df);
    CHECK(sol == doctest::Approx(stud).epsilon(1e-9));
  }
  
  TEST_CASE("bool sample_nonlinear_solver" * doctest::description("nonlinear solver")) {
    auto f = data.f; auto df = data.df;
    auto newt_scalar_step = [&f, &df] (double x) -> double {
      return mod_newt_step_scalar(x, f, df);
    };
    const double x_ex = std::tan(data.a);
    auto errf = [x_ex] (double& x) {
      return std::abs(x - x_ex);
    };
    
    double sol = 5, stud = 5;
    sample_nonlinear_solver(newt_scalar_step, sol, errf);
    sample_nonlinear_solver_TEST(newt_scalar_step, stud, errf);
    CHECK(sol == doctest::Approx(stud).epsilon(1e-9));
  }
  
  TEST_CASE("void mod_newt_ord" * doctest::description("scalar convergence order")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
  
  TEST_CASE("Vector mod_newt_step_system" * doctest::description("system modified Newton step")) {
    auto A = data.A; auto c = data.c;
    auto fn = [&A, &c](const Eigen::VectorXd& x){
      Eigen::VectorXd tmp =  A*x + c.cwiseProduct(x.array().exp().matrix()).eval();
      return tmp;
    };
    auto dfn = [&A, &c](const Eigen::VectorXd& x){
      Eigen::MatrixXd C = A;
      Eigen::VectorXd temp = c.cwiseProduct(x.array().exp().matrix());
      C += temp.asDiagonal();
      return C;
    };
    Eigen::VectorXd zeros = Eigen::VectorXd::Zero(5);
    const auto sol = mod_newt_step_system(zeros, fn, dfn);
    const auto stud = mod_newt_step_system_TEST(zeros, fn, dfn);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-9));
  }
  
  TEST_CASE("Eigen::VectorXd mod_newt_sys" * doctest::description("system iteration solver")) {
    const Eigen::VectorXd sol = mod_newt_sys(data.A, data.c);
    const Eigen::VectorXd stud = mod_newt_sys_TEST(data.A, data.c);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-9));
  }
  
  TEST_CASE("void mod_newt_sys_test" * doctest::description("system convergence order")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}


