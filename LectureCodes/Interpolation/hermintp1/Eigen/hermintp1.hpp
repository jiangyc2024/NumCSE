# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "hermloceval.hpp"

using Eigen::VectorXd;

// Compute and plot the cubic Hermite interpolant of the function \texttt{f} in the nodes \texttt{t}
// using weighted averages according to \eqref{pwintp:AverageSlopes} as local slopes 
template <class Function>
void hermintp1(Function& f, const VectorXd& t) {
  const unsigned n = t.size();
  // length of intervals between nodes: \Blue{$h_i := t_{i+1} - t_i$}
  const VectorXd h = t.tail(n - 1) - t.head(n - 1); 

  // evaluate f at the nodes and for plotting
  VectorXd y(n),
           tplot = VectorXd::LinSpaced(500, t(0), t(n - 1)),
           yplot(500);
  for (unsigned i = 0; i < n; ++i) y(i) = f(t(i));
  for (unsigned j = 0; j < 500; ++j) yplot(j) = f(tplot(j));

  // slopes of piecewise linear interpolant
  const VectorXd delta = (y.tail(n - 1) - y.head(n - 1)).cwiseQuotient(h);
  VectorXd c(n);
  c(0) = delta(0); c(n - 1) = delta(n - 2);
  // slopes from weighted averages, see \eqref{pwintp:AverageSlopes}
  for (unsigned i = 1; i < n - 1; ++i) {
    c(i) = ( h(i)*delta(i - 1) + h(i - 1)*delta(i) )/( t(i + 1) - t(i - 1) );
  }
  
  mgl::Figure fig;
  fig.title("Hermite interpolation");
  fig.legend(0, 0);
  fig.plot(t, y, " ko").label("Data points");
  fig.plot(tplot, yplot).label("f(x)");

  // compute and plot the Hermite interpolant with slopes \texttt{c}
  for (unsigned j = 0; j < n - 1; ++j) {
    VectorXd vx = VectorXd::LinSpaced(100, t(j), t(j + 1)),
             px;
    hermloceval(vx, t(j), t(j+1), y(j), y(j+1), c(j), c(j+1), px);
    fig.plot(vx, px, "r");
  }
  // manually adding label for pw hermite interpolants
  fig.addlabel("Piecw. cubic interpolation polynomial", "r");

  // plot segments indicating the slopes \Blue{$c_i$}
  for (unsigned k = 1; k < n - 1; ++k) {
    VectorXd t_sl(2), y_sl(2);
    t_sl << t(k) - 0.3*h(k - 1), t(k) + 0.3*h(k);
    y_sl << y(k) - 0.3*h(k - 1)*c(k), y(k) + 0.3*h(k)*c(k);
    fig.plot(t_sl, y_sl, "k");
  }
  // slope segments at beginning and end
  VectorXd t_sl_0(2), y_sl_0(2), 
           t_sl_n(2), y_sl_n(2);

  t_sl_0 << t(0), t(0) + 0.3*h(0);
  y_sl_0 << y(0), y(0) + 0.3*h(0)*c(0);
  t_sl_n << t(n - 1) - 0.3*h(n - 2), t(n - 1);
  y_sl_n << y(n - 1) - 0.3*h(n - 2)*c(n - 1), y(n - 1);
  fig.plot(t_sl_0, y_sl_0, "k");
  fig.plot(t_sl_n, y_sl_n, "k");

  fig.addlabel("Averaged slopes", "k");

  // save figure
  fig.save("hermintp1");
}
