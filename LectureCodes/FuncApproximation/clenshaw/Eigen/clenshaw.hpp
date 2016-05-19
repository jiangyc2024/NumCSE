# include <iostream>
# include <Eigen/Dense>
using Eigen::VectorXd;
using Eigen::MatrixXd;

// Clenshaw algorithm for evaluating \Blue{$p = \sum_{j=1}^{n+1}a_j T_{j-1}$} at points passed in vector \texttt{x}
// IN : \Blue{$\mathbf{a}$} in \Blue{$p = \sum_{j=1}^{n+1} a_j T_{j-1}
//      x = evaluation points
// OUT: values $p(x)$
VectorXd clenshaw(const VectorXd& a, const VectorXd& x) {
  const unsigned n = a.size() - 1; // degree of polynomial
  MatrixXd d(n + 1, x.size());
  for (unsigned c = 0; c < x.size(); ++c) {
    d.col(c) = a;
  }

  for (int j = n - 1; j > 0; --j) {
    d.row(j) += 2*x.transpose().cwiseProduct( d.row(j + 1) ); // see \eqref{eq:cstr}
    d.row(j - 1) -= d.row(j + 1);
  }

  return d.row(0) + x.transpose().cwiseProduct( d.row(1) );
}
