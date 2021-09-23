#ifndef PLOT_HPP
#define PLOT_HPP

#include <cassert>
#include <cmath>
#include <vector>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/* DO_NOT_CHANGE */
/**
 * @brief plot the vector sizes against the needed time to compute them
 * in loglog-plot
 */
void plot(std::vector<double> &vec_size, std::vector<double> &elap_time,
          std::vector<double> &elap_time_eff, const std::string &fig_name) {
  const unsigned int n = vec_size.size();
  // lines for the comparison of convergenz order
  std::vector<double> vec_size_lin(vec_size);
  std::vector<double> vec_size_pow3(vec_size);

  for (unsigned int i = 0; i < n; i++) {
    vec_size_lin[i] = vec_size[i] / vec_size[0] * elap_time[0];
    vec_size_pow3[i] =
        pow(vec_size[i], 3) / pow(vec_size[0], 3) * elap_time_eff[0];
  }

  // make sure the sizes match
  assert(vec_size.size() == elap_time.size() &&
         "vector sizes must be the same.");

  plt::figure();

  plt::loglog(vec_size, elap_time, "*", {{"label", "original"}});
  plt::loglog(vec_size, elap_time_eff, "*", {{"label", "effective"}});
  plt::loglog(vec_size, vec_size_lin, "k--", {{"label", "O(n)"}});
  plt::loglog(vec_size, vec_size_pow3, "--",
              {{"label", "O(n³)"}, {"color", "grey"}});
  plt::legend("best");

  plt::xlabel("Vector size (n)");
  plt::ylabel("Time [s]");
  plt::title("Comparison between timings");

  // note figname needs to have the right path: which is './cx_out/figname'
  plt::savefig(fig_name);
}

#endif
