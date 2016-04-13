# include <cmath>
# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include <timer.hpp>
# include "ANipoleval.hpp"
# include "ipolyeval.hpp"
# include "intpolyval.hpp"
# include "intpolyval_lag.hpp"

void print(const std::vector<double>& v) {
  for (const auto& val : v) {
    std::cout << val << " ";
  }
  std::cout << "\n";
}

int main() {
  // function to interpolate
  auto f = [](const Eigen::VectorXd& x){ return x.cwiseSqrt(); };

  const unsigned min_deg = 3, max_deg = 200;

  // checksum to prevent compiler from optimizing function calls away
  Eigen::VectorXd buffer;
  std::vector<double> t1, t2, t3, t4, N;

  // n = increasing polynomial degree
  for (unsigned n = min_deg; n <= max_deg; ++n) {

    double checksum;
    const Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(n, 1, n),
                          y = f(t);
    // ANipoleval takes a double as argument
    const double x = n*drand48(); // drand48 returns random double \Blue{$\in [0,1]$}
    // all other functions take a vector as argument
    const Eigen::VectorXd xv = n*Eigen::VectorXd::Random(1); 

    std::cout << "Degree = " << n << "\n";
//    timer<> aitken, ipol, intpol, intpol_lag;
    Timer aitken, ipol, intpol, intpol_lag;

    // do the same many times and choose the best result
    aitken.start();
    for (unsigned i = 0; i < 100; ++i)
      checksum += ANipoleval(t, y, x);                         aitken.lap();
    aitken.stop(); 
    t1.push_back(aitken.mean());
    std::cout << checksum << "\n";
    checksum = 0;

    ipol.start();
    for (unsigned i = 0; i < 100; ++i) 
      ipoleval(t, y, xv, buffer); checksum += buffer(0);       ipol.lap(); 
    t2.push_back(ipol.mean());
    std::cout << checksum << "\n";
    checksum = 0;

    intpol.start();
    for (unsigned i = 0; i < 100; ++i) 
      intpolyval(t, y, xv, buffer); checksum += buffer(0);     intpol.lap();
    t3.push_back(intpol.mean()); 
    std::cout << checksum << "\n";
    checksum = 0;

    intpol_lag.start();
    for (unsigned i = 0; i < 100; ++i) 
      intpolyval_lag(t, y, xv, buffer); checksum += buffer(0); intpol_lag.lap();
    t4.push_back(intpol_lag.mean());
    std::cout << checksum << "\n";
    checksum = 0;


    N.push_back(n);
    }

  print(t1);
  print(t2);
  print(t3);
  print(t4);

  mgl::Figure fig;
  fig.setlog(false, true);
  fig.plot(N, t1, "r").label("Aitken-Neville scheme");
  fig.plot(N, t2, "m").label("Naive polyfit + polyval");
  fig.plot(N, t3, "k").label("Barycentric formula");
  fig.plot(N, t4, "b").label("Lagrange polynomials");
  fig.xlabel("Polynomial degree");
  fig.ylabel("Computational time [s]");
  fig.legend(1,0);
  fig.save("pevaltime.eps");

  return 0;
}


      
      
