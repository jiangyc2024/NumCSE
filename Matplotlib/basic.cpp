#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#include <map>
#include <vector>

// Implements f(x) = x(x - 1) for x=-2...(n-2)
// Returns x, f(x)
std::pair<std::vector<double>, std::vector<double>> generate_data(const int n) {
  std::vector<double> x(n), y(n);
  for (int i = 0; i < n; ++i) {
    x[i] = i - 2;
    y[i] = x[i] * (x[i] - 1);
  }
  return std::make_pair(x, y);
}

int main() {
  // Get the data to plot
  const unsigned n = 20;
  auto data = generate_data(n);
  std::vector<double> x = data.first;
  std::vector<double> y = data.second;

  // Plot
  plt::figure();
  plt::title("Basic example");                          // plot title
  plt::plot(x, x);                                      // simple straight line
  plt::plot(x, y, {{"ls", "--"}, {"label", "x(x-1)"}}); // dashed with a label
  plt::grid(true);                                      // add a grid
  plt::xlabel("x");                                     // x axis label
  plt::ylabel("f(x)");                                  // y axis label
  plt::legend();                                        // activate legend
  plt::save("basic.pdf");

  // Log plot
  plt::figure();
  plt::title("Loglog plot: $f(x) = x(x-1)$");
  plt::loglog(x, y, "go-"); // python style definition of the linestyle
  plt::save("loglog.pdf");

  plt::show();

  return 0;
}
