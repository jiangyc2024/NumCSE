/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: December 2022
 * Homework problem "RK-SSMs and Discrete Evolutions", prb:rkssmeo
 */

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace ExplRKSSM {

// A class providing a general explicit Runge-Kutta single-step method definied
// through it Butcher tableaux.
/* SAM_LISTING_BEGIN_1 */
template <typename RHSFunctor>
class ExplRKSSMEvolOp {
 public:
  ExplRKSSMEvolOp(RHSFunctor& f, const Eigen::MatrixXd& A,
                  const Eigen::VectorXd& b)
      : f_(f), A_(A), b_(b) {}
  ExplRKSSMEvolOp(const ExplRKSSMEvolOp&) = default;
  ExplRKSSMEvolOp(ExplRKSSMEvolOp&&) = default;
  ExplRKSSMEvolOp& operator=(const ExplRKSSMEvolOp&) = default;
  ExplRKSSMEvolOp& operator=(ExplRKSSMEvolOp&&) = default;
  virtual ~ExplRKSSMEvolOp() = default;

  template <typename State>
  [[nodiscard]] State operator()(double h, const State& y) const;

 private:
  RHSFunctor& f_;      // Right-hand-side function $\Vf$
  Eigen::MatrixXd A_;  // Butcher matrix $\FA$
  Eigen::VectorXd b_;  // Weight vector $\Vb$
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <typename RHSFunctor>
template <typename State>
State ExplRKSSMEvolOp<RHSFunctor>::operator()(double h, const State& y) const {
  unsigned int s = A_.cols();
  assert(s == A_.rows());
  assert(s == b_.size());

  State y1{y};
  std::vector<State> k{s, y};
  for (int i = 0; i < (int)s; ++i) {
    for (int j = 0; j < i; ++j) {
      k[i] += h * A_(i, j) * k[j];
    }
    k[i] = f_(k[i]);
    y1 += h * b_[i] * k[i];
  }
  return y1;
}
/* SAM_LISTING_END_2 */

}  // namespace ExplRKSSM

int main(int /*argc*/, char** /*argv*/) {
  Eigen::MatrixXd A(4, 4);
  A.setZero();
  A(1, 0) = 0.5;
  A(2, 1) = 0.5;
  A(3, 2) = 1.0;
  Eigen::VectorXd b = Eigen::Vector4d(1 / 6.0, 2 / 6.0, 2 / 6.0, 1 / 6.0);

  auto f = [](Eigen::Vector2d& y) -> Eigen::Vector2d {
    return Eigen::Vector2d(y[0] * (1.0 - y[1]), y[1] * (y[0] - 1.0));
  };

  ExplRKSSM::ExplRKSSMEvolOp evolop(f, A, b);

  Eigen::Vector2d y1 = evolop(0.1, Eigen::Vector2d(2.0, 2.0));
  std::cout << "y1 = " << y1.transpose() << std::endl;
  return 0;
}
