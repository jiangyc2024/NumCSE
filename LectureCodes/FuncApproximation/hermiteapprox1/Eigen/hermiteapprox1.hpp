# include <iostream>
# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "hermloceval.hpp"
# include "polyfit.hpp"
# include "feval.hpp"
using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
void append(VectorXd&, const VectorXd&); // forward declaration of append()

// Investigation of interpolation error norms for cubic Hermite interpolation of \Blue{$f$} (handle \texttt{f})
// on \Blue{$[a,b]$} with linearly averaged slopes according to \eqref{pwintp:AverageSlopes}.
// \texttt{N} gives the maximum number of mesh intervals
template <class Function, class Derivative>
void hermiteapprox(const Function& f, const Derivative& df, const double a, const double b, const unsigned N) {
  std::vector<double> l2err, linferr, h; // save error and stepwidths in these vectors 
  for (unsigned j = 2; j <= N; ++j) {
    // \texttŧ{xx} is the fine mesh on which the error norms are computed
    VectorXd xx(1); xx << a;     
    // \texttt{val} contains the hermite approximated values in \texttt{xx}
    VectorXd val(1); val << f(a);

    VectorXd t = VectorXd::LinSpaced(j, a, b); // mesh nodes
    VectorXd y = feval(f, t); // feval returns f evaluated at t
    VectorXd c = feval(df, t); // coefficients for Hermit polynomial representation

    for (unsigned k = 0; k < j - 1; ++k) {
      // local evaluation nodes in interval \Blue{$[t_k, t_{k+1}]$}
      VectorXd vx = VectorXd::LinSpaced(100, t(k), t(k+1)),
               locval;
      // evaluate the hermite interpolant at the \texttt{vx}, using \cref{hermloceval}
      hermloceval(vx, t(k), t(k+1), y(k), y(k+1), c(k), c(k+1), locval);
      // do not append the first value as that one is already in the vector
      append(xx, vx.tail(99));
      append(val, locval.tail(99));
    }

    // difference between exact function and interpolant
    VectorXd d = feval(f, xx) - val; 
    const double L2 = d.lpNorm<2>(), // \Blue{$L^2$} norm
                 Linf = d.lpNorm<Eigen::Infinity>(); // \Blue{$L^{\infty}$} norm
    l2err.push_back(L2); // save L2 error
    linferr.push_back(Linf); // save Linf error
    h.push_back( (b - a)/j ); // save current meshwidth
  }

  // compute estimates for algebraic orderns of convergence 
  // using linear regression on half the data points

  // L2Log contains the logarithmus of the last M l2err values,
  // LinfLog and hLog the other according values
  const unsigned M = unsigned( N / 2 );
  VectorXd L2Log(M), LinfLog(M), hLog(M);
  for (unsigned n = 0; n < M; ++n) {
    // *(h.end() - n - 1) = h[end - n]
    L2Log(M - n - 1) = std::log( *(l2err.end() - n - 1) );
    LinfLog(M - n - 1) = std::log( *(linferr.end() - n - 1) );
    hLog(M - n - 1) = std::log( *(h.end() - n - 1) );
  }

  // use linear regression
  VectorXd pI = polyfit( hLog, LinfLog, 1 );
  std::cout << "Algebraic convergence rate of Infinity norm: " << pI(0) << "\n";
  VectorXd pL2 = polyfit( hLog, L2Log, 1 );
  std::cout << "Algebraic convergence rate of L2 norm: " << pL2(0) << "\n";

  mgl::Figure fig;
  fig.setlog( true, true ); // double logarithmic plot
  fig.title( "Hermite approximation" );
  fig.xlabel( "Meshwidth h" );
  fig.ylabel( "Norm of interpolation error" );
  fig.legend();
  fig.plot( h, l2err, " +r" ).label( "L^2 Error" );
  fig.plot( h, linferr, " +b" ).label( "L^{\\infty} Error" );

  // plot linear regression lines
  fig.plot( std::vector<double>({h[0], h[N-2]}), std::vector<double>({std::pow(h[0], pL2(0))*std::exp(pL2(1)), std::pow(h[N-2], pL2(0))*std::exp(pL2(1))}), "m" );
  fig.plot( std::vector<double>({h[0], h[N-2]}), std::vector<double>({std::pow(h[0], pI(0))*std::exp(pI(1)), std::pow(h[N-2], pI(0))*std::exp(pI(1))}), "k" );
  fig.save( "hermiperravgsl" );
}

// Appends a Eigen::VectorXd to another Eigen::VectorXd
void append(VectorXd& x, const VectorXd& y) {
  x.conservativeResize(x.size() + y.size());
  x.tail(y.size()) = y;
}
/* SAM_LISTING_END_0 */
