#include <vector>

#include "pchi.hpp"

using namespace Eigen;

int main() {
  std::vector<double> N_nodes, h, err_reconstr, err_zero;

  for (unsigned int i = 4; i <= 512; i = i << 1) {
    N_nodes.push_back(i);
    h.push_back(10. / i);
  }
  err_zero = fppchipConvergence();
  err_reconstr = rspchipConververgence();

  for (int i = 0; i < h.size(); ++i) {
    std::cout << N_nodes[i] << " " << err_zero[i] << " " << err_reconstr[i]
              << std::endl;
  }

  // Error plot
  plt::figure();
  plt::title("Error VS no. of nodes");
  plt::xlabel("h");
  plt::ylabel("max |f(t) - s(t)|");
  plt::loglog(h, err_reconstr, "b", {{"label", "$s_{reconstr}$"}});
  plt::loglog(h, err_zero, "r", {{"label", "$s_{zero}$"}});

  std::vector<double> vec_size_lin(h.size());
  std::vector<double> vec_size_pow3(h.size());

  for (unsigned int i = 0; i < h.size(); i++) {
    vec_size_lin[i] = h[i] / h[0];
    vec_size_pow3[i] = pow(h[i], 3) / pow(h[0], 3);
  }

  plt::loglog(h, vec_size_lin, "k--", {{"label", "$O(n)$"}});
  plt::loglog(h, vec_size_pow3, "k--", {{"label", "$O(n^3)$"}});

  plt::legend();
  plt::savefig("./cx_out/pchi_conv.png");
}
