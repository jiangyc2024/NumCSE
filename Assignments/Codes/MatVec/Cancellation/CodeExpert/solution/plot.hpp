
#ifndef PLOT_HPP
#define PLOT_HPP

#include <cmath>
#include <Eigen/Dense>
#include <vector>

# include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace Eigen;


void py_plot( ArrayXd& h, ArrayXd& ex, ArrayXd& g1, ArrayXd& g2 )
{
  ArrayXd err2 = (g2-ex).abs();
  ArrayXd err1 = (g1-ex).abs();

  
  plt::figure();
  plt::loglog(h, h, "--", {{"label", "O(h)"}});
  plt::loglog(h, err2, "--", {{"label", "g2-ex"}});
  plt::loglog(h, err1, "--", {{"label", "g1-ex"}});

  plt::title("Error of approximation of f'(x_0)");
  plt::xlabel("h");
  plt::ylabel("| f'(x_0) - g_i(x_0, h) |");
  plt::legend("best");
  plt::savefig("./cx_out/error_cancellation.png");


};
    
#endif 
