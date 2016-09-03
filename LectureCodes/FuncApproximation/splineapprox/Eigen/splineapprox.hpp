# include <iostream>
# include <cmath> // needed for std::log()
# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "polyfit.hpp" // provides polyfit(), see NumCSE/Utils
# include "feval.hpp" // provides feval(), see NumCSE/Utils
# include "completespline.hpp" // provides spline()
using Eigen::VectorXd;

// Study of interpolation error norms for \emph{complete} cubic spline interpolation of \Blue{$f$}
// on equidistant knot set.
template <class Function, class Derivative>
void splineapprox(Function& f, Derivative& df, const double a, const double b, const unsigned N, const std::string plotname) {
  const unsigned K = int( (b-a)/0.00025 ); // no. of evaluation points for norm evaluation
  const VectorXd x = VectorXd::LinSpaced(K, a, b), // mesh for norm evaluation
                fv = feval(f, x);

  std::vector<double> errL2, errInf, h;
  VectorXd v, // save result of spline approximation here
           t; // spline knots
  for (unsigned j = 3; j <= N; ++j) {
    t = VectorXd::LinSpaced(j, a, b); // spline knots
    VectorXd ft = feval(f, t);

    // compute complete spline imposing exact first derivative at the endpoints 
    // spline: takes set of interpolation data (t, y) and slopes at endpoints,
    //         and returns values at points x
    v = spline(t, ft, df(a), df(b), x);
    // compute error norms 
    h.push_back( (b - a)/j ); // save current mesh width
    errL2.push_back( (fv - v).lpNorm<2>() );
    errInf.push_back( (fv - v).lpNorm<Eigen::Infinity>() );
  }

  // compute algebraic orders of convergence using polynomial fit
  const unsigned n = h.size();
  VectorXd hLog(n), errL2Log(n), errInfLog(n);
  for (unsigned i = 0; i < n; ++i) {
    hLog(i) = std::log(h[i]);
    errL2Log(i) = std::log(errL2[i]);
    errInfLog(i) = std::log(errInf[i]);
  }

  VectorXd p = polyfit(hLog, errL2Log, 1);
  std::cout << "L2 norm convergence rate: " << p(0) << "\n";
  p = polyfit(hLog, errInfLog, 1);
  std::cout << "Supremum norm convergence rate: " << p(0) << "\n";

  // plot interpolation
  mgl::Figure fig;
  fig.title("Spline interpolation " + plotname);
  fig.plot(t, feval(f, t), " m*").label("Data points");
  fig.plot(x, fv, "b").label("f");
  fig.plot(x, v, "r").label("Cubic spline interpolant");
  fig.legend(1,0);
  fig.save("interp_" + plotname);

  // plot error
  mgl::Figure err;
  err.title("Spline approximation error " + plotname);
  err.setlog(true, true);
  err.plot(h, errL2, "r;").label("L^2 norm");
  err.plot(h, errInf, "b;").label("L^\\infty norm");
  fig.legend(1, 0);
  err.save("approx_" + plotname);
}
