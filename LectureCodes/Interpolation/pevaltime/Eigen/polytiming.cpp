# include <cmath>
# include <vector>
# include <Eigen/Dense>
# include <figure.hpp>
# include "timer.h"
# include "ANipoleval.hpp"
# include "ipolyeval.hpp"
# include "intpolyval.hpp"
# include "intpolyval_lag.hpp"

int main() {
  // function to interpolate
  auto f = [](const Eigen::VectorXd& x){ return x.cwiseSqrt(); };

  const unsigned min_deg = 3, max_deg = 200;

  Eigen::VectorXd buffer;
  std::vector<double> t1, t2, t3, t4, N;
  
  // Number of repeats for each eval
  const int repeats = 100;

  // n = increasing polynomial degree
  for (unsigned n = min_deg; n <= max_deg; n++) {

    const Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(n, 1, n),
                          y = f(t);

    // ANipoleval takes a double as argument
    const double x = n*drand48(); // drand48 returns random double \Blue{$\in [0,1]$}
    // all other functions take a vector as argument
    const Eigen::VectorXd xv = n*Eigen::VectorXd::Random(1); 

    std::cout << "Degree = " << n << "\n";
    Timer aitken, ipol, intpol, intpol_lag;

    // do the same many times and choose the best result
    aitken.start();
    for (unsigned i = 0; i < repeats; ++i){
      ANipoleval(t, y, x); aitken.lap(); }
    t1.push_back(aitken.min());

    ipol.start();
    for (unsigned i = 0; i < repeats; ++i) {
      ipoleval(t, y, xv, buffer); ipol.lap(); }
    t2.push_back(ipol.min());

    intpol.start();
    for (unsigned i = 0; i < repeats; ++i) {
      intpolyval(t, y, xv, buffer); intpol.lap(); }
    t3.push_back(intpol.min()); 

    intpol_lag.start();
    for (unsigned i = 0; i < repeats; ++i) {
      intpolyval_lag(t, y, xv, buffer); intpol_lag.lap(); }
    t4.push_back(intpol_lag.min());

    N.push_back(n);
    }

  mgl::Figure fig;
  fig.setlog(false, true);
  fig.title("Timing for single-point evaluation");
  fig.plot(N, t1, "r").label("Aitken-Neville scheme");
  fig.plot(N, t2, "m").label("Polyfit + Polyval");
  fig.plot(N, t3, "k").label("Barycentric formula");
  fig.plot(N, t4, "b").label("Lagrange polynomials");
  fig.xlabel("Polynomial degree");
  fig.ylabel("Computational time [s]");
  fig.legend(0,1);
  fig.save("pevaltime.eps");

  return 0;
}


      
      
