# include <Eigen/Dense>
# include <complex>

namespace intpolyval {


using Eigen::VectorXcd;

inline
// IN:  \texttt{t}: vector of nodes \Blue{$t_0, \ldots, t_n$}
//      \texttt{y}: vector of data \Blue{$y_0, \ldots, y_n$}
//      \texttt{x}: vector of evaluation points \Blue{$x_1, \ldots, x_N$}
// OUT: \texttt{p}: interpolant evaluated at x 
void intpolyval(const VectorXcd& t, const VectorXcd& y, const VectorXcd& x, VectorXcd& p) {
  const unsigned n = t.size(); // no. of interpolation nodes = deg. of polynomial $-1$
  const unsigned N = x.size(); // no. of evaluation points

  p = VectorXcd(N); // resizing 

  // Precompute the weights \Blue{$\lambda_i$} with effort \Blue{$O(n^2)$}
  VectorXcd lambda(n);
  for (unsigned k = 0; k < n; ++k) {
    // little workaround: cannot subtract a vector from a scalar 
    // -> multiply scalar by vector of ones
    lambda(k) = 1./( (t(k)*VectorXcd::Ones(k) - t.head(k)).prod() * 
                     (t(k)*VectorXcd::Ones(n - k - 1) - t.tail(n - k - 1)).prod() );
  }
  // Compute quotient of weighted sums  of \Blue{$\frac{\lambda_i}{t - t_i}$}, effort \Blue{$O(n)$}
  for (unsigned i = 0; i < N; ++i) {
    VectorXcd z = (x(i)*VectorXcd::Ones(n) - t);

    // check if we want to evaluate at a node <-> avoid division by zero
    auto ptr = std::find(z.begin(), z.end(), std::complex<double>(0,0));
    if (ptr != z.end()) { // if ptr = z.end no zero was found
      p(i) = y(ptr - z.begin()); // ptr - z.begin gives the position of the zero
    }
    else {
      const VectorXcd mu = lambda.cwiseQuotient(z);
      p(i) = (mu.cwiseProduct(y)).sum()/mu.sum();
    }
  }
}


} //namespace intpolyval