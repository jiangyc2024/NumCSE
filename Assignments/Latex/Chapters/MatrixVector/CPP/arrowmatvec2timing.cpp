# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "arrowmatvec.hpp"
# include "arrowmatvec2.hpp"
# include "timer.h"
# include "timer1.h"


/**************************************************************************/
# include <iostream>

void libtimer() {
  std::cout << "timer.h:\n";
  for (unsigned n = 4; n <= 4000; n *= 2){
    // create test input 
    Eigen::VectorXd a = Eigen::VectorXd::Random(n),
                    d = Eigen::VectorXd::Random(n),
                    x = Eigen::VectorXd::Random(n),
                    y1, y2;
    timer<> t1;
    t1.start(); arrowmatvec(d, a, x, y1); t1.stop();
    std::cout << "n = " << n << ": " << t1.elapsed().count()/1.e9 << "  ";

    timer<> t2;
    t2.start(); arrowmatvec2(d, a, x, y2); t2.stop();
    std::cout << t2.elapsed().count()/1.e9 <<  "\n";
  }
}

void libtimer1() {
  std::cout << "timer1.h:\n";
  for (unsigned n = 4; n <= 4000; n *= 2){
    // create test input 
    Eigen::VectorXd a = Eigen::VectorXd::Random(n),
                    d = Eigen::VectorXd::Random(n),
                    x = Eigen::VectorXd::Random(n),
                    y1, y2;
    Timer t1;
    t1.start(); arrowmatvec(d, a, x, y1); t1.stop();
    std::cout << "n = " << n << ": " << t1.duration() << "  ";
    
    Timer t2;
    t2.start(); arrowmatvec2(d, a, x, y2); t2.stop();
    std::cout << t2.duration() <<  "\n";
  }
}

void sanity_check() {
  // loops which are also used to plot in the "timing()" funtion (below)
  libtimer();
  libtimer1();  

  libtimer1();
  libtimer();
}
/**************************************************************************/

void timing() {
  std::vector<double> evals, res1, res2; // save timings and vector sizes here
  for (unsigned n = 4; n <= 8192; n *= 2){
    // create test input 
    Eigen::VectorXd a = Eigen::VectorXd::Random(n),
                    d = Eigen::VectorXd::Random(n),
                    x = Eigen::VectorXd::Random(n),
                    y1, y2;
    timer<> t1, t2;
    t1.start(); arrowmatvec(d, a, x, y1); t1.stop();
    t2.start(); arrowmatvec2(d, a, x, y2); t2.stop();
    evals.push_back(n); // save vector sizes
    res1.push_back(t1.elapsed().count()/1e9); // time gets measured in nanoseconds
    res2.push_back(t2.elapsed().count()/1e9); 
  }

  mgl::Figure fig;
  fig.title("Timings arrowmatvec2");
  //fig.ranges(2, 9000, 1e-6, 1e3);
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
  sanity_check();
  //timing();
  return 0;
}
