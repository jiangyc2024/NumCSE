
#ifndef PLOT_HPP
#define PLOT_HPP

#include <Eigen/Dense>
#include <cmath>
#include <vector>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/* DO_NOT_CHANGE */
/* @brief plot the vector sizes against the needed time to compute them
 * in loglog-plot
 */
void plot(std::vector<double> &vec_size, std::vector<double> &elap_time1,
          std::vector<double> &elap_time2, const std::string &fig_name,
          const std::string &label1, const std::string &label2) {
  unsigned int n = vec_size.size();
  // lines for the comparison of convergenz order
  std::vector<double> vec_size_lin(vec_size);
  std::vector<double> vec_size_pow2(vec_size);

  for (unsigned int i = 0; i < n; i++) {
    vec_size_lin[i] = vec_size[i] / vec_size[0] * elap_time2[0] / 100.0;
    vec_size_pow2[i] = 
	pow(vec_size[i], 2) / pow(vec_size[0], 2) * elap_time2[0] / 100.0;
  }

  // make sure the sizes match
  assert(vec_size.size() == elap_time1.size() &&
         "vector sizes must be the same.");

  plt::figure();
  // Multiplication with ten to shift the plots
  plt::loglog(vec_size, elap_time1, "+", {{"label", label1}});
  plt::loglog(vec_size, elap_time2, "+", {{"label", label2}});
  plt::loglog(vec_size, vec_size_lin, "k--",
              {{"label", "O(n)"}, {"color", "grey"}});
  plt::loglog(vec_size, vec_size_pow2, "-",
              {{"label", "0(n²)"}, {"color", "black"}});
  plt::legend("best");

  plt::xlabel("Vector size (n)");
  plt::ylabel("Time [s]");
  plt::title("Comparison between timings");

  // note figname needs to have the right path: which is './cx_out/figname'
  plt::savefig(fig_name);
}

#endif
