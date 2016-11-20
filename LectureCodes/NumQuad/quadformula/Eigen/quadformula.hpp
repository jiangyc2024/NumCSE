// Generic numerical quadrature routine implementing \eqref{eq:quadform}: 
// \texttt{f} is a handle to a function, e.g. as lambda function
// \texttt{c}, \texttt{w} pass quadrature nodes \Blue{$\qn_j\in [a,b]$}, and weights \Blue{$\qw_j>0$}
// in a Eigen::VectorXd
template <class Function>
double quadformula(Function& f, const Eigen::VectorXd& c, const Eigen::VectorXd& w) {
  const std::size_t n = c.size();
  double I = 0;
  for (std::size_t i = 0; i < n; ++i) I += w(i)*f(c(i));
  return I;
}
