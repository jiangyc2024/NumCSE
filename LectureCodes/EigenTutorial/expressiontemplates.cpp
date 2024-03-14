// **********************************************************************
// Eigen tutorial codes: http://eigen.tuxfamily.org/dox/GettingStarted.html
// **********************************************************************

#include "timetable.hpp"
#include <Eigen/Dense>
using Eigen::MatrixXd;

/*

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

void init(size_t s) {
  M = MatrixXd::Random(s, s);
  A = MatrixXd::Random(s, s);
  B = MatrixXd::Random(s, s);
}

*/

int main() {

  MatrixXd M;
  MatrixXd A;
  MatrixXd B;
  MatrixXd R;

  // Invert D as expression template
  // Not forcing the type of certain eigen expressions enables internal
  // optimizations in eigen with expression templates (by lazy evaluation)
  auto fastInvertDiagonal = [&]() {
    auto D = M.diagonal().asDiagonal();
    R = D.inverse();
  };

  // Invert D (forced) as dense matrix
  auto slowInvertDiagonal = [&]() {
    MatrixXd D = M.diagonal().asDiagonal();
    R = D.inverse();
  };

  // Add three matrices with intermediate evaluations
  auto slowAddThree = [&]() { 
    R = 3 * M;
    R += 4 * A;
    R += 5 * B; 
  };

  // Add three matrices as one line
  // Large expressions give eigen more opportunities for optimization, in this
  // case eigen can compile the whole computation into one loop instead of 3 (or
  // more)
  auto fastAddThree = [&]() { R = 3 * M + 4 * A + 5 * B; };

  auto init = [&](size_t s) {
    M = MatrixXd::Random(s, s);
    A = MatrixXd::Random(s, s);
    B = MatrixXd::Random(s, s);
  };

  const std::vector<size_t> sizes = {10, 20, 40, 80, 160, 320, 640};

  timeTable(
      sizes,
      {slowInvertDiagonal, fastInvertDiagonal, slowAddThree, fastAddThree},
      init, {"size", "slow inv.", "fast inv.", "slow add", "fast add"});
}