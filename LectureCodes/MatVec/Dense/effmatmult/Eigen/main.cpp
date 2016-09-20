///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>

#include <Eigen/Dense>
#include <figure/figure.hpp>

#include "effmatmult.hpp"


int main(){
  Eigen::MatrixXd timings = dottenstiming();
  std::cout << timings << std::endl;
  //Plotting
  mgl::Figure fig;
  fig.setFontSize(5);
  fig.title("Timings for rank 1 matrix-vector multiplications");
  fig.setlog(true, true);
  fig.plot(timings.col(0),timings.col(1), " +b").label("efficient evaluation");
  fig.plot(timings.col(0),timings.col(2)," ^r").label("slow evaluation");
  fig.plot(timings.col(0),timings.col(0).array().pow(1).matrix() * timings(timings.rows()-1,1) / (std::pow(timings(timings.rows()-1,0),1)), "l;").label("O(n)");
  fig.plot(timings.col(0),timings.col(0).array().pow(2).matrix() * timings(timings.rows()-1,2) / (std::pow(timings(timings.rows()-1,0),2)), "h;").label("O(n^2)");
  fig.xlabel("vector length n");
  fig.ylabel("time [s]");
  fig.legend(0.05,0.95);
  fig.save("dottenstiming");
  return 0;
}
