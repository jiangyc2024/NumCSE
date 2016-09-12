#include <iostream>
// INTERNAL ONLY
#include <chrono>
// INTERNAL ONLY END

#include <Eigen/Dense>

using namespace Eigen;

//! @brief Build an "arrow matrix"
//! Given vectors $a$ and $b$, returns A*A*x in $y$, where A is build from a,d
//! @param[in] d A n-dimensional vector
//! @param[in] a A n-dimensional vector
//! @param[in] x A n-dimensional vector
//! @param[out] y The vector y = A*A*x
void efficient_arrow_matrix_2_times_x(const VectorXd & d, const VectorXd & a,
                              const VectorXd & x, VectorXd & y) {
  assert(d.size() == a.size() && a.size() == x.size() &&
         "Vector size must be the same!");
  int n = d.size();

  // SOLUTION ONLY

  // Notice that we can compute (A*A)*x more efficiently using
  // A*(A*x). This is in fact performing two matrix vetor multiplications
  // instead of a more expensive matrix-matrix multiplication.
  // Therefore, as first step, we need a way to efficiently compute A*x

  // This function computes A*x
  // Notice that A = D + H + V, s.t. A*x = D*x + H*x + V*x
  // D*x can be rewritten as d*x componentwise
  // H*x is zero, except at the last component
  // V*x
  auto A_times_x = [] (const VectorXd & x) {
    // This takes care of the diagonal (D*x)
    VectorXd Ax = ( d.array() * x.array() ).matrix();

    // Now H*x only affects the last component of A*x
    Ax(n-1) += a.head(n - 1) * x.head(n-1);

    // Now V*x is equal to the vector (a(0)*x(n-1), ..., a(n-2)*x(n-1), 0)
    Ax.head(n-1) +=  x(n-1) * a.head(n-1);

    return Ax;
  }

  y = A_times_x(A_times_x(x));

  // END SOLUTION ONLY
}

//! @brief Build an "arrow matrix"
//! Given vectors $a$ and $b$, returns A*A*x in $y$, where A is build from a,d
//! @param[in] d A n-dimensional vector
//! @param[in] a A n-dimensional vector
//! @param[in] x A n-dimensional vector
//! @param[out] y The vector y = A*A*x
void arrow_matrix_2_times_x(const VectorXd & d, const VectorXd & a,
                              const VectorXd & x, VectorXd & y) {
  assert(d.size() == a.size() && a.size() == x.size() &&
         "Vector size must be the same!");
  int n = d.size();

  // BEGIN SOLUTION COMMENT
  // In this lines, we extract the blocks used to construct the matrix A
  // END SOLUTION COMMENT
  VectorXd d_head = d.head(n-1);
  VectorXd a_head = a.head(n-1);
  MatrixXd d_diag = d_head.asDiagonal();

  MatrixXd A(n,n);

  // BEGIN SOLUTION COMMENT
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
  // END SOLUTION COMMENT
  A << d_diag,             a_head,
       a_head.transpose(), d(n-1);

  y = A*A*x;
}

// INTERNAL ONLY
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

  for(int n = 2; n < 2048; n = n << 1) {
    VectorXd d = VectorXd::Random(n);
    VectorXd a = VectorXd::Random(n);
    VectorXd x = VectorXd::Random(n);
    VectorXd y(n);
    duration_t elapsed = timing([&a, &d, &x, &y] () { arrow_matrix_2_times_x(d,a,x,y); }, 1);
    std::cout << elapsed.count() << " ns" << std::endl;
  }

}
// END INTERNAL ONLY

int main(void) {
  // VectorXd a = {1,2,3,4,5};
  // VectorXd d = {1,3,4,5,6};
  // VectorXd x = {-5,4,6,-8,5};

  time_arrow_matrix();
}
