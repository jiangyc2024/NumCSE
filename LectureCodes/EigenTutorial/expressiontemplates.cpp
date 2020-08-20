// **********************************************************************
// Eigen tutorial codes: http://eigen.tuxfamily.org/dox/GettingStarted.html
// **********************************************************************

#include <Eigen/Dense>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;
using namespace std::chrono;

MatrixXd M, A, B, R;

// Helper function to time an action over N repetitions
template <class Action>
size_t time(Action&& a, size_t const N) {
  size_t minimum = ~0;
  for (int i = 0; i < N; ++i) {
    auto start = high_resolution_clock::now();
    a();
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    minimum = std::min(minimum, (size_t)duration.count());
  }

  return minimum;
}

// Invert D as expression template
// Not forcing the type of certain eigen expressions enables internal
// optimizations in eigen with expression templates (by lazy evaluation)
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
// Large expressions give eigen more opportunities for optimization, in this
// case eigen can compile the whole computation into one loop instead of 3 (or
// more)
void fastAddThree() { R = 3 * M + 4 * A + 5 * B; }

// Add three matrices with intermediate evaluations
void slowAddThree() {
  R = 3 * M;
  R += 4 * A;
  R += 5 * B;
}

int main() {
  vector<double> sizes = {10, 20, 40, 80, 160, 320, 640};
  size_t N(5);

  vector<double> slowInv(7), fastInv(7), slowAdd(7), fastAdd(7);
  cout << setw(5) << "size" << setw(10) << "slow inv." << setw(10)
       << "fast inv." << setw(10) << "slow add" << setw(10) << "fast add"
       << endl;

  for (size_t i = 0; i < sizes.size(); ++i) {
    auto s = (size_t)sizes[i];
    M = MatrixXd::Random(s, s);
    A = MatrixXd::Random(s, s);
    B = MatrixXd::Random(s, s);

    slowInv[i] = time(&slowInvertDiagonal, N);
    fastInv[i] = time(&fastInvertDiagonal, N);
    slowAdd[i] = time(&slowAddThree, N);
    fastAdd[i] = time(&fastAddThree, N);

    cout << setw(5) << s << setw(10) << slowInv[i] << setw(10) << fastInv[i]
         << setw(10) << slowAdd[i] << setw(10) << fastAdd[i] << endl;
  }
}