# include <iostream>
# include <complex>
# include <cmath> // M\_PI
# include <limits> // numeric\_limits<double>::max()
# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "timer.h"
# include "trigipequid.hpp"
# include "trigpolycoeff.hpp"
using Eigen::VectorXd;

// Runtime comparison between efficient ($\to$ Code~\ref{trigipequid}) and direct computation
// ($\to$ Code~\ref{trigpolycoeff} of coefficients of trigonoetric interpolation polynomial in
// equidistant points.
void trigipequidtiming() {
  const int Nruns = 5;
  std::array< std::vector<double>, 3 > times; // times = [n, t1, t2]

  for (int n = 10; n <= 500; n += 10) {
    std::cout << "n = " << n << "\n";
    const int N = 2*n + 1;
    const VectorXd t = VectorXd::LinSpaced(N, 0, 1 - 1./N),
          y = (2*M_PI*t).array().cos().exp().matrix();

    Timer t1, t2;
    VectorXd a, b;
    VectorXcd ac, bc;
    std::complex<double> checksum; // to avoid optimizing function calls away
    for (int k = 0; k < Nruns; ++k) {
      t1.start(); trigpolycoeff(t, y, a, b); t1.lap();
      checksum += a(0) + b(0);
      t2.start(); trigipequid(y, ac, bc); t2.lap();
      checksum += ac(0) + bc(0);
    }

    times[0].push_back(n); // save number of evaluations
    times[1].push_back(t1.min()); // save time of trigpolycoeff
    times[2].push_back(t2.min()); // save time of trigipequid
  }

  mgl::Figure fig;
  fig.setlog(true, true);
  fig.plot(times[0], times[1], " b+").label("trigpolycoeff");
  fig.plot(times[0], times[2], " r*").label("trigipequid");
  fig.xlabel("n");
  fig.ylabel("runtime [s]");
  fig.legend(0, 1);
  fig.save("trigipequidtiming");
}

int main() {
  trigipequidtiming();
  return 0;
}
      
