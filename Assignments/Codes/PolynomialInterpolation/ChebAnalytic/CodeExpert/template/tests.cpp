#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"
#include <random>

// includes for test data
#include <Eigen/Dense>

TEST_SUITE("ChebAnalytic") {

  TEST_CASE("double lengthEllipse" *
            doctest::description("Approximate curve length of an ellipse")) {
    unsigned int n = 100;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.001, 20.0);
    double rho = dis(gen);
    double result = lengthEllipse_TEST(rho, n);
    // calculate length of ellipse using series expansion
    std::vector<long double> coeff{1.L/4, 1.L/64, 1.L/256, 25.L/16384, 49.L/65536, 441.L/1048576};
    long double ref = 1.0;
    long double h = 1.0L/(rho*rho*rho*rho);
    long double hp = h;
    for (Eigen::Index i = 0; i < 6; i++) {
      ref += coeff[i]*hp;
      hp *= h;
    }
    ref *= M_PI*rho;
    CHECK(abs(result - ref) == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("std::pair<double, double> bestBound" *
            doctest::description("calculate the approximation of $M = \\min_\\rho"
                                 " \\frac{\\vert\\gamma_\\rho\\vert \\max_{z \\in \\gamma_p}"
                                 " \\vert f(z) \\vert}{\\pi d([-1,1], \\gamma)}$")) {
    auto sol_M_rho = bestBound();
    double sol_M = sol_M_rho.first;
    double sol_rho = sol_M_rho.second;
    auto stud_M_rho = bestBound_TEST();
    double stud_M = stud_M_rho.first;
    double stud_rho = stud_M_rho.second;

    CHECK(abs(sol_M - stud_M) == doctest::Approx(0.).epsilon(1e-6));
    CHECK(abs(sol_rho - stud_rho) == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("std::pair<std::vector<double>, std::vector<double>> compareErrorAndBound" *
            doctest::description("Tabulate the estimated upperbound and actual error for Chebychev interpolation")) {
    auto n_max = 20;
    auto sol_err_up = compareErrorAndBound(n_max);
    std::vector<double> sol_err_v = sol_err_up.first;
    std::vector<double> sol_upbd_v = sol_err_up.second;
    auto stud_err_up = compareErrorAndBound_TEST(n_max);
    std::vector<double> stud_err_v = stud_err_up.first;
    std::vector<double> stud_upbd_v = stud_err_up.second;
    
    REQUIRE(sol_err_v.size() == stud_err_v.size());
    REQUIRE(sol_upbd_v.size() == stud_upbd_v.size());
    Eigen::Map<Eigen::VectorXd> sol_err(sol_err_v.data(),sol_err_v.size());
    Eigen::Map<Eigen::VectorXd> stud_err(stud_err_v.data(),stud_err_v.size());
    Eigen::Map<Eigen::VectorXd> sol_upbd(sol_upbd_v.data(),sol_upbd_v.size());
    Eigen::Map<Eigen::VectorXd> stud_upbd(stud_upbd_v.data(),stud_upbd_v.size());
    CHECK((sol_err-stud_err).norm() == doctest::Approx(0.).epsilon(1e-9));

    CHECK((sol_upbd-stud_upbd).norm() == doctest::Approx(0.).epsilon(1e-9));
  }

}

