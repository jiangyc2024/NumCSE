#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

#include "timer.h"
#if INTERNAL
#include <chrono>
#include <figure/figure.hpp>
#endif // INTERNAL

using namespace Eigen;

/* @brief Build an "arrow matrix"
 * Given vectors $a$ and $b$, returns A*A*x in $y$, where A is build from a,d
 * @param[in] d A n-dimensional vector
 * @param[in] a A n-dimensional vector
 * @param[in] x A n-dimensional vector
 * @param[out] y The vector y = A*A*x
 */
/* SAM_LISTING_BEGIN_1 */
void efficient_arrow_matrix_2_times_x(const VectorXd & d,
                                      const VectorXd & a,
                                      const VectorXd & x,
                                      VectorXd & y) {
  assert(d.size() == a.size() && a.size() == x.size() &&
         "Vector size must be the same!");
  int n = d.size();

  #if SOLUTION
  // Notice that we can compute (A*A)*x more efficiently using
  // A*(A*x). This is in fact performing two matrix vetor
  // multiplications
  // instead of a more expensive matrix-matrix multiplication.
  // Therefore, as first step, we need a way to efficiently
  // compute A*x

  // This function computes A*x, you can use it
  // by calling A\_times\_x(x);
  // This is the syntax for lambda functions: notice the extra
  // [variables] code. Each variable written within [] brakets
  // will be caputred (i.e. seen) inside the lambda function
  // Without \&, a copy of the variable will be performed

  // Notice that A = D + H + V, s.t. A*x = D*x + H*x + V*x
  // D*x can be rewritten as d*x componentwise
  // H*x is zero, except at the last component
  // V*x is similar is only affected by the last component of x
  auto A_times_x = [&a, &d, n] (const VectorXd & x) {
    // This takes care of the diagonal (D*x)
    // Notice: we use d.array() to tell Eigen to treat
    // a vector an an array. As a result: each operation
    // is performed componentwise.
    VectorXd Ax = ( d.array() * x.array() ).matrix();

    // H*x only affects the last component of A*x
    // This is a dot product between a and x with the last
    // component removed
    Ax(n-1) += a.head(n - 1).dot(x.head(n-1));

    // V*x is equal to the vector
    // (a(0)*x(n-1), ..., a(n-2)*x(n-1), 0)
    Ax.head(n-1) +=  x(n-1) * a.head(n-1);

    return Ax;
  };

  // <=> y = A*A*x
  y = A_times_x(A_times_x(x));

  #endif // SOLUTION
}
/* SAM_LISTING_END_1 */

/* @brief Build an "arrow matrix"
 * Given vectors $a$ and $b$, returns A*A*x in $y$, where A is build from a,d
 * @param[in] d A n-dimensional vector
 * @param[in] a A n-dimensional vector
 * @param[in] x A n-dimensional vector
 * @param[out] y The vector y = A*A*x
 */
/* SAM_LISTING_BEGIN_0 */
void arrow_matrix_2_times_x(const VectorXd & d, const VectorXd & a,
                              const VectorXd & x, VectorXd & y) {
  assert(d.size() == a.size() && a.size() == x.size() &&
         "Vector size must be the same!");
  int n = d.size();

  #if SOLUTION
  // In this lines, we extract the blocks used to construct the matrix A.
  #endif // SOLUTION
  VectorXd d_head = d.head(n-1);
  VectorXd a_head = a.head(n-1);
  MatrixXd d_diag = d_head.asDiagonal();

  MatrixXd A(n,n);

  #if SOLUTION
  // We build the matrix A using the "comma initialization": each expression separated
  // by a comma is a "block" of the matrix we are building.
  // d\_diag is the top left (n-1)x(n-1) block
  // a\_head is the top right vertical vector
  // a\_head.transpose() is the bottom left horizontal vector
  // d(n-1) is a single element (a 1x1 matrix), on the bottom right corner
  // This is how the matrix looks like:
  // A = | D   | a      |
  //     |-----+--------|
  //     | a\^T | d(n-1) |
  #endif // SOLUTION
  A << d_diag,             a_head,
       a_head.transpose(), d(n-1);

  y = A*A*x;
}
/* SAM_LISTING_END_0 */

#if INTERNAL
using duration_t = std::chrono::nanoseconds;

template <class Function>
duration_t timing(const Function & F, int repeats = 10) {

  using time_point_t = std::chrono::high_resolution_clock::time_point;

  duration_t min_elapsed;
  for(int r = 0; r < repeats; r++) {

    time_point_t start = std::chrono::high_resolution_clock::now();

    F();

    duration_t elapsed = std::chrono::duration_cast<duration_t>(
                     std::chrono::high_resolution_clock::now() - start);
    min_elapsed = r == 0 ? elapsed : std::min(elapsed, min_elapsed);
  }

  return min_elapsed;

}

void time_arrow_matrix() {

  std::vector<double> sizes, timings, timings_efficient;

  for(int n = (1 << 5); n < (1 << 12); n = n << 1) {
    VectorXd d = VectorXd::Random(n);
    VectorXd a = VectorXd::Random(n);
    VectorXd x = VectorXd::Random(n);
    VectorXd y(n);
    duration_t elapsed = timing([&a, &d, &x, &y] () {
        arrow_matrix_2_times_x(d,a,x,y);
      }, 1);
    duration_t elapsed_efficient = timing([&a, &d, &x, &y] () {
        efficient_arrow_matrix_2_times_x(d,a,x,y);
      }, 1);

    std::cout << n
              << std::scientific << std::setprecision(3)
              << std::setw(15) << elapsed.count()
              << std::setw(15) << elapsed_efficient.count()
              << std::endl;

    timings.push_back(elapsed.count() * 10e-9);
    timings_efficient.push_back(elapsed_efficient.count() * 10e-9);
    sizes.push_back(n);
  }

  mgl::Figure fig;
  fig.title("Timings of arrow\\_matrix\\_2\\_times\\_x");
  fig.ranges(2, 9000, 1e-8, 1e3);
  fig.setlog(true, true); // set loglog scale
  fig.plot(sizes, timings, " r+").label("runtime");
  fig.fplot("1e-9*x^3", "k|").label("O(n^3)");
  fig.xlabel("Vector size (n)");
  fig.ylabel("Time [s]");
  fig.legend(0, 1);
  fig.save("arrowmatvec_timing.eps");
  fig.save("arrowmatvec_timing.png");

  mgl::Figure fig2;
  fig2.title("Comparison of timings");
  fig2.ranges(2, 9000, 1e-8, 1e3);
  fig2.setlog(true, true); // set loglog scale
  fig2.plot(sizes, timings, " r+").label("original");
  fig2.plot(sizes, timings_efficient, " r+").label("efficient");
  fig2.fplot("1e-9*x^3", "k|").label("O(n^3)");
  fig2.fplot("1e-7*x", "k").label("O(n)");
  fig2.xlabel("Vector size (n)");
  fig2.ylabel("Time [s]");
  fig2.legend(0, 1);
  fig2.save("arrowmatvec_comparison.eps");
  fig2.save("arrowmatvec_comparison.png");
}
#endif // INTERNAL

void runtime_arrow_matrix() {
#if SOLUTION
  // sizes will contain the size of the matrix
  // timings will contain the runtimes in seconds
  /* SAM_LISTING_BEGIN_3 */
  std::vector<double> sizes, timings;

  std::cout << "n" << "\t" << "time" << std::endl;
  for (unsigned n = 4; n <= 2048; n *= 2){
    // Create test input using random vectors
    Eigen::VectorXd a = Eigen::VectorXd::Random(n),
                    d = Eigen::VectorXd::Random(n),
                    x = Eigen::VectorXd::Random(n),
                    y;
    Timer t;
    t.start();
    // Perform a single test
    // TODO: more tests
    arrow_matrix_2_times_x(d, a, x, y);
    efficient_arrow_matrix_2_times_x(d, a, x, y);
    t.stop();
    sizes.push_back(n); // save vector sizes
    timings.push_back(t.duration());

    std::cout << n << "\t" << t.duration() << std::endl;
  }
  /* SAM_LISTING_END_3 */
#else
  // TODO: your code here
#endif // SOLUTION
}

#if SOLUTION
void runtime_arrow_matrix_with_chrono() {
  std::vector<double> sizes, timings;
  std::cout << "n" << "\t" << "time" << std::endl;
  for (unsigned n = 4; n <= 2048; n *= 2){
    // Create test input using random vectors
    Eigen::VectorXd a = Eigen::VectorXd::Random(n),
                    d = Eigen::VectorXd::Random(n),
                    x = Eigen::VectorXd::Random(n),
                    y;
    Timer t;
    t.start();
    // Perform a single test
    // TODO: more tests
    arrow_matrix_2_times_x(d, a, x, y);
    efficient_arrow_matrix_2_times_x(d, a, x, y);
    t.stop();
    sizes.push_back(n); // save vector sizes
    timings.push_back(t.duration());
    std::cout << n << "\t" << t.duration() << std::endl;
  }
}
#endif // SOLUTION

int main(void) {
  VectorXd a(5);
  a << 1., 2., 3., 4., 5.;
  VectorXd d(5);
  d <<1., 3., 4., 5., 6.;
  VectorXd x(5);
  x << -5., 4., 6., -8., 5.;
  VectorXd yi, ye;

  arrow_matrix_2_times_x(a,d,x,yi);
  efficient_arrow_matrix_2_times_x(a,d,x,ye);

  double err = (yi - ye).norm();

  std::cout << "Error: " << err << std::endl;

#if INTERNAL
  time_arrow_matrix();
#endif // INTERNAL

  runtime_arrow_matrix();
#if SOLUTION
  runtime_arrow_matrix_with_chrono();
#endif // SOLUTION

  double eps = std::numeric_limits<double>::denorm_min();
  exit(err < eps);
}
