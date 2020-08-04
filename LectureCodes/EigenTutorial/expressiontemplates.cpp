// **********************************************************************
// Eigen tutorial codes: http://eigen.tuxfamily.org/dox/GettingStarted.html
// **********************************************************************

#include <Eigen/Dense>
#include <chrono>
#include <iostream>

using namespace Eigen;
using namespace std;
using namespace std::chrono;

MatrixXd M, A, B, R;

// Helper function to time an action over N repetitions
template <class Action>
void time(Action&& a, size_t const N) {
  auto start = high_resolution_clock::now();

  for (int i = 0; i < N; ++i) a();

  auto end = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  cout << "took " << duration.count() << " ms" << endl;
}

// Invert D as expression template
void fastInvertDiagonal() {
  auto D = M.diagonal().asDiagonal();
  R = D.inverse();
}

// Invert D (forced) as dense matrix
void slowInvertDiagonal() {
  MatrixXd D = M.diagonal().asDiagonal();
  R = D.inverse();
}

// Add three matrices as one line
void fastAddThree() { R = 3 * M + 4 * A + 5 * B; }

// Add three matrices with intermediate evaluations
void slowAddThree() {
  R = 3 * M;
  R += 4 * A;
  R += 5 * B;
}

int main() {
  size_t S(10);
  size_t N(1000);
  M = MatrixXd::Random(S, S);
  A = MatrixXd::Random(S, S);
  B = MatrixXd::Random(S, S);

  // Not forcing the type of certain eigen expressions enables internal
  // optimizations in eigen with expression templates (by lazy evaluation)
  cout << "fast invert diagonal ";
  time(&fastInvertDiagonal, N);
  cout << "slow invert diagonal ";
  time(&slowInvertDiagonal, N);

  // Large expressions give eigen more opportunities for optimization, in this
  // case eigen can compile the whole computation into one loop instead of 3 (or
  // more)
  cout << "fast add three matrices ";
  time(&fastAddThree, N);
  cout << "slow add three matrices ";
  time(&slowAddThree, N);
}