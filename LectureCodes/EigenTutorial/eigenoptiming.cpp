/* Demonstraction code for course Numerical Methods  for CSE, ETH Zurich
   Tim9ing of linear algebra operations in Eigen
   @author Ralf Hiptmair
   @date August 2020
*/

#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std::chrono;

/* SAM_LISTING_BEGIN_T */
void eigenoptiming(std::vector<int> &&nvals, const char *filename,
                   unsigned int Nrep = 10) {
  assert(!nvals.empty());
  // Vector for collecting runtimes
  std::vector<std::tuple<int, double, double, double>> times{};
  for (int n : nvals) {
    // Allocated test vectors an matrices
    const Eigen::VectorXd v = Eigen::VectorXd::Random(n);
    const Eigen::MatrixXd M = Eigen::MatrixXd::Random(n, n);
    const Eigen::MatrixXd X = Eigen::MatrixXd::Random(n, n);
    // Result variables
    Eigen::VectorXd y = Eigen::VectorXd::Zero(n);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    double s = 0.0;
    // Vriables for recording mimimal runtimes
    double t_dotp = 1.0E9;
    double t_mv = 1.0E9;
    double t_mm = 1.9E9;
    for (unsigned int j = 0; j < Nrep; ++j) {
      {
        auto start = high_resolution_clock::now();
        s += v.dot(v);
        auto end = high_resolution_clock::now();
        auto duration = (duration_cast<microseconds>(end - start)).count();
        t_dotp = (duration < t_dotp) ? duration : t_dotp;
      }
      {
        auto start = high_resolution_clock::now();
        y += M * v;
        auto end = high_resolution_clock::now();
        auto duration = (duration_cast<microseconds>(end - start)).count();
        t_mv = (duration < t_mv) ? duration : t_mv;
      }
      {
        auto start = high_resolution_clock::now();
        A += M * X;
        auto end = high_resolution_clock::now();
        auto duration = (duration_cast<microseconds>(end - start)).count();
        t_mm = (duration < t_mm) ? duration : t_mm;
      }
    }
    times.push_back({n, t_dotp, t_mv, t_mm});
  }
  // Output to file
  std::ofstream outfile(filename);
  for (const auto &data : times) {
    auto [n, t_dotp, t_mv, t_mm] = data;
    outfile << n << " , " << t_dotp << " , " << t_mv << " , " << t_mm
            << std::endl;
  }
  outfile.close();
}

/* SAM_LISTING_END_T */

int main(int /*argc*/, char ** /*argv*/) {
  eigenoptiming({10, 20, 40, 80, 160, 320, 640, 1280}, "eigenoptimings.csv");
  return 0;
}
