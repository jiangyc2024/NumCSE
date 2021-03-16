#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

int main(int argc, char **argv) {
  /* SAM_LISTING_BEGIN_1 */
  int n = 10;

#if SOLUTION
  std::vector<double> e(n), loge(n);

  e[0] = 1;
  loge[0] = std::log(e[0]);
  e[1] = 0.8;
  loge[1] = std::log(e[1]);
  for (int k = 1; k < n; ++k) {
    e[k + 1] = e[k] * std::sqrt(e[k - 1]);
    loge[k + 1] = std::log(e[k + 1]);

    std::cout << (loge[k + 1] - loge[k]) / (loge[k] - loge[k - 1]) << std::endl;
  }
#else
  // TODO: approximate order
#endif
  /* SAM_LISTING_END_1 */
}
