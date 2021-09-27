#include <Eigen/Dense>
#include <fstream>
#include <vector>

#include "radioactive.hpp"

int main() {
  std::ifstream data_file("decay.txt");
  if (!data_file.is_open()) {
    std::cerr << "Cannot open decay.txt\n";
    return 1;
  }

  std::vector<double> t_data, m_data;
  for (double t_i, m_i; data_file >> t_i >> m_i;) {
    t_data.push_back(t_i);
    m_data.push_back(m_i);
  }
  Eigen::Map<Eigen::ArrayXd> t(t_data.data(), t_data.size()),
      m(m_data.data(), m_data.size());

  // Implements the function object for F s.t. F(x) = 0 from subtask b
  /* SAM_LISTING_BEGIN_0 */
  auto F = [&t, &m](const Eigen::Array4d& x) -> Eigen::VectorXd {
    Eigen::VectorXd f = Eigen::VectorXd::Zero(t.size());
    // TODO: (9-11.d) Use Eigen's array operations which use componentwise
    // arithmetic to compute F
    // START
    const double a0 = x[0], b0 = x[1], l1 = x[2], l2 = x[3];
    f = (-l2 * t).exp() * b0 +
        (l1 / (l2 - l1)) * ((-l1 * t).exp() - (-l2 * t).exp()) * a0 - m;
    // END
    return f;
  };
  /* SAM_LISTING_END_0 */

  // Implements the Jacobian of F from subtask c
  /* SAM_LISTING_BEGIN_1 */
  auto DF = [&t](const Eigen::Array4d& x) -> Eigen::MatrixX4d {
    Eigen::MatrixX4d df = Eigen::MatrixX4d::Zero(t.size(), 4);
    // TODO: (9-11.d) Use Eigen's array operations which use componentwise
    // arithmetic to compute the Jacobian of F
    // START
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
    // END
    return df;
  };
  /* SAM_LISTING_END_1 */

  Eigen::Vector4d x(1., 1., 1., 0.1);
  std::vector<double> gn_update = GaussNewton(x, F, DF);

  plot(t, m, x, gn_update);
}
