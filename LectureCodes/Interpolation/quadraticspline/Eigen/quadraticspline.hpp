
# include <Eigen/Dense>

using Eigen::VectorXd;
// choice of extra knots \Blue{$p_i \in \left]t_{i-1}, t_i\right[$}
// for given nodeset \texttt{t} and coefficients \texttt{c}
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd extra_knots(const Eigen::VectorXd& t, const Eigen::VectorXd &y,
			    const Eigen::VectorXd& c) {
  const unsigned n = t.size();
  assert((n == y.size()) && (n == c.size()));
  Eigen::VectorXd p = (t[n - 1] + 1.0)*Eigen::VectorXd::Ones(n - 1);
  for (unsigned j = 0; j < n - 1; ++j) {
    // Safe comparison, just to avoid division by zero exception
    // If both slopes a close $\cob{p_j}$ will certainly lie outside the node interval 
    if (c[j] != c[j + 1]) {
      p[j] = (y[j+1] - y[j] + t[j]*c[j] - t[j+1]*c[j+1]) / (c[j+1]- c[j]);
    }
    if (p[j] < t[j] || p[j] > t[j+1]) {
      p[j] = 0.5*(t[j] + t[j+1]);
    }
  }
  return p;
}

/* SAM_LISTING_END_1 */
