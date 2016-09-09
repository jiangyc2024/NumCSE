# include <cmath>
# include <vector>
# include <Eigen/Dense>

using Eigen::VectorXd;

// Computation of \Hyperlink{LEBESGUE}{Lebesgue constant} of polynomial interpolation 
// with knots \Blue{$t_i$} passed in the vector \texttt{t} based on \eqref{eq:IPN1}. 
// \texttt{N} specifies the number of sampling points for the approximate 
// computation of the maximum norm of the Lagrange polynimial on the interval \Blue{$[-1,1]$}.
double lebesgue(const VectorXd& t, const unsigned& N) {
  const unsigned n = t.size();

  // compute denominators of normalized Lagrange polynomials relative to the nodes t
  VectorXd den(n);
  for (unsigned i = 0; i < n; ++i) {
    VectorXd tmp(n - 1);
    // Note: comma initializer can't be used with vectors of length 0
    if (i == 0) tmp = t.tail(n - 1);
    else if (i == n - 1) tmp = t.head(n - 1);
    else tmp << t.head(i), t.tail(n - (i + 1));
    den(i) = (t(i) - tmp.array()).prod();
  }

  double l = 0; // return value
  for (unsigned j = 0; j < N; ++j) {
    const double x = -1 + j*(2./N); // sampling point for \Blue{$\NxLinf{\cdot}{[-1,1]}$}
    double s = 0;
    for (unsigned k = 0; k < n; ++k) {
      // \texttt{v} provides value of normalized Lagrange polynomials
      VectorXd tmp(n - 1);
      if (k == 0) tmp = t.tail(n - 1);
      else if (k == n - 1) tmp = t.head(n - 1);
      else tmp << t.head(k), t.tail(n - (k + 1));
      double v = (x - tmp.array()).prod()/den(k);
      s += std::abs(v); // sum over modulus of the polynomials
    }
    l = std::max(l, s); // maximum of sampled values
  }
  return l;
}
