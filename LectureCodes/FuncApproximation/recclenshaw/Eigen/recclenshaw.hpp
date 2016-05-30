# include <Eigen/Dense>
using Eigen::VectorXd;

// Recursive evaluation of a polynomial \Blue{$p = \sum_{j=1}^{n+1}a_j T_{j-1}$} at point \texttt{x}, see \eqref{eq:cstr}
// IN : Vector of coefficients a
//      evaluation point x
// OUT: Value at point x
double recclenshaw(const VectorXd& a, const double x) {
  const int n = a.size() - 1;
  std::cout << n << "\n";
  double y;
  if (n < 2) {
    // \Blue{$y = a_2*x + a_1$}
    y = x*a(1) + a(0); 
  }
  else {
    VectorXd new_a(n);
    new_a << a.head(n - 2), a(n - 2) - a(n), a(n - 1)+ 2*x*a(n);
    y = recclenshaw(new_a, x); // recurse
  }
  return y;
}
