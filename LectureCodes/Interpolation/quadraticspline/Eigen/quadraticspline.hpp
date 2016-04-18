# include <Eigen/Dense>
using Eigen::VectorXd;

// choice of extra knots \Blue{$p_i \in \left]t_{i-1}, t_i\right[
// for given nodeset \texttt{t} and coefficients \texttt{c}
void extra_knots(const VectorXd& t, const VectorXd& c, VectorXd& p) {
  const unsigned n = t.size();
  p = (t(n - 1) - 1)*VectorXd::Ones(n - 1);
  for (unsigned j = 0; j < n - 1; ++j) {
    if (c(j) != c(j + 1)) {
      p(j) = (y(j+1) - y(j) + t(j)*c(j) - t(j+1)*c(j+1)) / (c(j+1)- c(j));
    }
    if (p(j) < t(j) || p(j) > t(j+1)) {
      p(j) = 0.5*(t(j) + t(j+1));
    }
  }
}
