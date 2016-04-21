# include <Eigen/Dense>

using mat_t = Eigen::MatrixXd;
using vec_t = Eigen::VectorXd;

// returns the values of the first n - 1 legendre polynomials in x
// as columns of the matrix L 
void legendre(const unsigned& n, const vec_t& x, mat_t& L) {
  L = mat_t::Ones(n,n); // \Blue{$p_0(x) = 1$}
  L.col(1) = x; // \Blue{$p_1(x) = x$}
  for (unsigned j = 1; j < n - 1; ++j) {
    // \Blue{$p_{j+1}(x) = \frac{2j + 1}{j + 1} x p_{j}(x) - \frac{j}{j + 1} p_{j - 1}(x)$} \cref{eq:Legpol}
    V.col(j + 1) = (2.*j + 1)/(j + 1.)*V.col(j - 1).cwiseProduct(x) - j/(j + 1.)*V.col(j - 1);
  }
}
