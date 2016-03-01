template <class Function>
double quadformula(Function& f, const Eigen::VectorXd& c, const Eigen::VectorXd& w) {
// Generic numerical quadrature routine implementing \eqref{eq:quadform}: 
// \texttt{f} is a handle to a function, e.g. as lambda function
// \texttt{c}, \texttt{w} pass quadrature nodes \Blue{$\qn_j\in [a,b]$}, and weights \Blue{$\qw_j>0$}
// in a Eigen::VectorXd

  const long n = c.size();
  double I = 0;
  for (long i = 0; i < n; ++i) {
    I += w(i)*f(c(i));
  }

  return I;
}
