# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "arrowmatvec.hpp"
# include "arrowmatvec2.hpp"
# include "timer1.h"


void timing() {
  std::vector<double> evals, res1, res2; // save timings and vector sizes here
  for (unsigned n = 4; n <= 8192; n *= 2){
    // create test input 
    Eigen::VectorXd a = Eigen::VectorXd::Random(n),
                    d = Eigen::VectorXd::Random(n),
                    x = Eigen::VectorXd::Random(n),
                    y1, y2;
    Timer t1, t2;
    t1.start(); arrowmatvec(d, a, x, y1); t1.stop();
    t2.start(); arrowmatvec2(d, a, x, y2); t2.stop();
    evals.push_back(n); // save vector sizes
    res1.push_back(t1.duration()); // time gets measured in nanoseconds
    res2.push_back(t2.duration()); 
  }

  mgl::Figure fig;
  fig.title("Timings arrowmatvec2");
  fig.ranges(2, 9000, 1e-8, 1e3);
  fig.setlog(true, true); // set loglog scale
  fig.plot(evals, res1, " r+").label("arrowmatvec");
  fig.plot(evals, res2, " g+").label("arrowmatvec2");
  fig.fplot("1e-9*x^3", "k|").label("O(n^3)"); 
  fig.fplot("1e-7*x", "k").label("O(n)");
  fig.xlabel("Vector size n");
  fig.ylabel("Time [s]");
  fig.legend(0, 1);
  fig.save("arrowmatvec2timing.eps");
  fig.save("arrowmatvec2timing.png");
}

int main() {
  timing();
  return 0;
}
