# include <Eigen/Dense>

using Eigen::VectorXd;

// using forward declaration of the function pchipend, implementation below
double pchipend(const double&, const double&, const double&, const double&);

void pchipslopes(const VectorXd& t, const VectorXd& y, VectorXd& c) {
  // Calculation of local slopes \Blue{$c_i$} for shape preserving cubic Hermite interpolation, see \eqref{mteq:lim}, \eqref{mteq:hm}
  // \texttt{t}, \texttt{y} are vectors passing the data points
  const unsigned n = t.size();
  const VectorXd h = t.tail(n - 1) - t.head(n - 1),
                 delta = (y.tail(n - 1) - y.head(n - 1)).cwiseQuotient(h); // linear slopes
  c = VectorXd::Zero(n);
  
  // compute reconstruction slope according to \eqref{mteq:hm}
  for (unsigned i = 0; i < n - 2; ++i) {
    if (delta(i)*delta(i + 1) > 0) {
      const double w1 = 2*h(i + 1) + h(i),
                   w2 = h(i + 1) + 2*h(i);
      c(i + 1) = (w1 + w2)/(w1/delta(i) + w2/delta(i + 1));
    }
  }
  // Special slopes at endpoints, beyond \eqref{mteq:hm}
  c(0) = pchipend(h(0), h(1), delta(0), delta(1));
  c(n - 1) = pchipend(h(n - 2), h(n - 3), delta(n - 2), delta(n - 3));
}

double pchipend(const double& h1, const double& h2, const double& del1, const double& del2) {
  // Non-centered, shape-preserving, three-point formula
  double d = ((2*h1 + h2)*del1 - h1*del2)/(h1 + h2);
  if (d*del1 < 0) {
    d = 0;
  }
  else if (del1*del2 < 0 && std::abs(d) > std::abs(3*del1)) {
    d = 3*del1;
  }
  return d;
}
