# ifndef POLYFIT_HPP
# define POLYFIT_HPP

# include <Eigen/Dense>
# include <Eigen/QR>

Eigen::VectorXd polyfit(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const unsigned& order) {
  // A = [1 x_1 x_1^2 ... ]
  //     [ ...        ... ]
  //     [1 x_n x_n^2 ... ]
  Eigen::MatrixXd A = Eigen::MatrixXd::Ones(x.size(), order + 1);
  for (unsigned j = 1; j < order + 1; ++j) {
    A.col(j) = A.col(j - 1).cwiseProduct(x);
  }
  Eigen::VectorXd coeffs = A.householderQr().solve(y);

  return coeffs.reverse();
}

# endif
