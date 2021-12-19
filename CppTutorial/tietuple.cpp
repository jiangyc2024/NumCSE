/***********************************************************************
 *                                                                     *
 * Demo code                                                           *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: R.H.                                                        *
 * Date: May 2018                                                      *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

#include <cmath>
// Header for basic IO
#include <iostream>
// Provides random acccess container class
#include <vector>
// Provides algorithms operating on generic containers
#include <algorithm>
// Miscellaneous utilities
#include <tuple>

using std::cout;
using std::endl;

template <typename T>
std::tuple<T, T, std::vector<T>> extcumsum(const std::vector<T> &v) {
  // Local summation variable captured by reference by lambda function
  T sum{};
  // temporary vector for returning cumulative sum
  std::vector<T> w{};
  // cumulative summation
  std::transform(v.cbegin(), v.cend(), back_inserter(w), [&sum](T x) {
    sum += x;
    return (sum);
  });
  return (std::tuple<T, T, std::vector<T>>(
      *std::min_element(v.cbegin(), v.cend()),
      *std::max_element(v.cbegin(), v.cend()), w));
}

void tietuple() {
  // initialize a vector from an initializer list
  std::vector<double> v({1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8});
  // Variables for return values
  double minv = NAN;
  double maxv = NAN;      // Extremal elements
  std::vector<double> cs; // Cumulative sums
  std::tie(minv, maxv, cs) = extcumsum(v);
  cout << "min = " << minv << ", max = " << maxv << endl;
  cout << "cs = [ ";
  for (double x : cs) {
    cout << x << ' ';
  }
  cout << "]" << endl;
}

void compoundreturn() {
  // initialize a vector from an initializer list
  std::vector<double> v({1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8});
  // Definition of variables and assignment of return values all at once
  auto [minv, maxv, cs] = extcumsum(v);
  cout << "min = " << minv << ", max = " << maxv << endl;
  cout << "cs = [ ";
  for (double x : cs) {
    cout << x << ' ';
  }
  cout << "]" << endl;
}

int main(int /*argc*/, char ** /*argv*/) {
  tietuple();
  compoundreturn();
  return 0;
}
