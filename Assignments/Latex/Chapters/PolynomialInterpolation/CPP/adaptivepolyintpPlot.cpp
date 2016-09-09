# include <cmath>
# include <figure/figure.hpp>
# include "./adaptivepolyintp.hpp"

int main() {
  // declare test functions 
  auto f1 = [](double t) { return std::sin(std::exp(2*t)); };
  auto f2 = [](double t) { return std::sqrt(t)/(1 + 16*t*t); };

  // test interval 
  const double a = 0, b = 1;

  // get interpolation nodes and print runtimes
  const unsigned N = 1000; // no. of sampling points
  const double tol = 1e-6; // tolerance
  Eigen::VectorXd tf1, tf2, // nodes
                  ef1, ef2; // errors

  adaptivepolyintp(f1, a, b, tol, N, tf1, ef1);
  adaptivepolyintp(f2, a, b, tol, N, tf2, ef2);

  // plot
  mgl::Figure fig;
  fig.title("Errors at each step");
  fig.setlog(false, true);
  fig.xlabel("No. of interpolation nodes");
  fig.ylabel("max |f(t) - I_Tf(t)|");
  fig.plot(ef1, "r").label("f_1(t) = sin(e^{2t})");
  fig.plot(ef2, "b").label("f_2(t) = \\sqrt{t}/(1 + 16t^2)");
  fig.legend();
  fig.save("cvg");

  return 0;
}
