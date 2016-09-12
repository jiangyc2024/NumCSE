# include <string>
# include <cmath>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "trapezoidal.hpp"

void cvg(const double l, const double r, const bool logx, const std::string title, const std::string saveas) {

  // repeating the experiment for different parameters 'a'
  Eigen::VectorXd a(4); a << 0.5, 0.9, 0.95, 0.99;
  // applying trapezoidal rule for 1..N intervals
  const unsigned N = 30;

  mgl::Figure fig;
  fig.title(title);
  fig.xlabel("no. of quadrature nodes");
  fig.ylabel("|quadrature error|");
  if (logx) {
    fig.setlog(true, true); // loglog scale
  } else {
    fig.setlog(false, true); // semilogy scale
  }
  std::vector<std::string> plotstyles = {"b+", "r+", "g+", "c+"};
  for (unsigned i = 0; i < a.size(); ++i) {
    // we must define the function here, because the trapezoidal function requires a function
    // that takes only one argument! so we can't give the parameter 'a' as argument and 
    // therefore must define the function again for every different 'a' we use
    auto issin = [a, i](double x) {
      // f(x) = 1/sqrt(1 - a*sin(2*pi*x + 1)
      return 1/std::sqrt(1 - a(i)*std::sin(2*M_PI*x + 1));
    };

    std::vector<double> res; // saving results here
    res.reserve(N);
    // compute approximation for exact value by choosing a large no. of steps
    const double ex = trapezoidal(issin, l, r, N*200);

    for (unsigned n = 1; n <= N; ++n) {
      // computing error for 1..N steps
      res.push_back( std::abs(ex - trapezoidal(issin, l, r, n)) );
    }
    // plot error
    fig.plot(res, plotstyles[i]).label("a = " + std::to_string(a(i)));
  }

  fig.legend(0, 0); // activate legend and put it in bottom left corner
  fig.save(saveas);  // save plot

  return;
}

int main () {
  // test for non periodic (using loglog scale)
  cvg(0, 0.5, true, "Non periodic: on [0, 0.5]", "traperr2.eps");
  // test for periodic interval (using semilogy scale)
  cvg(0, 1, false, "Periodic: on [0,1]", "traperr1.eps");

  return 0;
}





