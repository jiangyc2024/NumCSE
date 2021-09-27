#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <vector>

#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("Radioactive") {
  TEST_CASE("std::vector<double> GaussNewton" *
            doctest::description("Gauss-Newton method")) {
    std::ifstream data_file("decay.txt");
    if (!data_file.is_open()) {
      std::cerr << "Cannot open decay.txt\n";
      REQUIRE(false);
    }

    std::vector<double> t_data, m_data;
    for (double t_i, m_i; data_file >> t_i >> m_i;) {
      t_data.push_back(t_i);
      m_data.push_back(m_i);
    }
    Eigen::Map<Eigen::ArrayXd> t(t_data.data(), t_data.size()),
        m(m_data.data(), m_data.size());

    auto F = [&t, &m](const Eigen::Array4d& x) -> Eigen::VectorXd {
      Eigen::VectorXd f = Eigen::VectorXd::Zero(t.size());
      const double a0 = x[0], b0 = x[1], l1 = x[2], l2 = x[3];
      f = (-l2 * t).exp() * b0 +
          (l1 / (l2 - l1)) * ((-l1 * t).exp() - (-l2 * t).exp()) * a0 - m;
      return f;
    };
    auto DF = [&t](const Eigen::Array4d& x) -> Eigen::MatrixX4d {
      Eigen::MatrixX4d df = Eigen::MatrixX4d::Zero(t.size(), 4);
      const double a0 = x[0], b0 = x[1], l1 = x[2], l2 = x[3];
      const double inv_d = 1. / (l2 - l1);
      Eigen::ArrayXd expl1 = (-l1 * t).exp();
      Eigen::ArrayXd expl2 = (-l2 * t).exp();
      Eigen::ArrayXd expd = expl1 - expl2;

      Eigen::VectorXd col1 = l1 * inv_d * expd;
      Eigen::VectorXd col2 = expl2;
      Eigen::VectorXd col3 =
          l2 * inv_d * inv_d * expd * a0 - l1 * t * inv_d * expl1 * a0;
      Eigen::VectorXd col4 = -t * expl2 * b0 - l2 * inv_d * inv_d * expd * a0 +
                             l1 * t * inv_d * expl2 * a0;
      df << col1, col2, col3, col4;
      return df;
    };

    Eigen::Vector4d x(1., 1., 1., 0.1);
    Eigen::Vector4d x_stud = x;

    // We only test if the correct result is produced.

    GaussNewton(x, F, DF, 1e-14);
    GaussNewton_TEST(x_stud, F, DF, 1e-14);

    CHECK((x - x_stud).norm() == doctest::Approx(0.).epsilon(1e-10));
  }

  TEST_CASE("void plot" * doctest::description("Visualization")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
