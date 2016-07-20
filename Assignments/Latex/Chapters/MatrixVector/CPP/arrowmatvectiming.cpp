# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "arrowmatvec.hpp"
# include "timer.h"

void timing() {
  std::vector<double> evals, res; // save timings and vector sizes here
  for (unsigned n = 4; n <= 8192; n *= 2){
    // create test input 
    Eigen::VectorXd a = Eigen::VectorXd::Random(n),
                    d = Eigen::VectorXd::Random(n),
                    x = Eigen::VectorXd::Random(n),
                    y;
    timer<> t;
    t.start(); arrowmatvec(d, a, x, y); t.stop();
    evals.push_back(n); // save vector sizes
    res.push_back(t.elapsed().count()/1e9); // time gets measured in nanoseconds
  }

  mgl::Figure fig;
  fig.title("Timings arrowmatvec");
  fig.ranges(2, 9000, 1e-6, 1e3);
  fig.setlog(true, true); // set loglog scale
  fig.plot(evals, res, " r+").label("runtime");
  fig.fplot("1e-9*x^3", "k|").label("O(n^3)"); 
  fig.xlabel("Vector size n");
  fig.ylabel("Time [s]");
  fig.legend(0, 1);
  fig.save("arrowmatvectiming.eps");
  fig.save("arrowmatvectiming.png");
}

int main() {
  timing();
  return 0;
}
