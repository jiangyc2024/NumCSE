# include <Eigen/Dense>
# include "meshgrid.hpp" // provided by NumCSE/Utils
using Eigen::MatrixXd; using Eigen::VectorXd;

void psf(const long L, MatrixXd& S) {
  VectorXd x = VectorXd::LinSpaced(2*L+1, -L, L);
  MatrixXd X,Y;
  meshgrid(x, x, X, Y);
  MatrixXd E = MatrixXd::Ones(2*L+1, 2*L+1);
  S = E.cwiseQuotient(E + X.cwiseProduct(X) + Y.cwiseProduct(Y));
  S /= S.sum();
}
