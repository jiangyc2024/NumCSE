// **********************************************************************
// Eigen tutorial codes: http://eigen.tuxfamily.org/dox/GettingStarted.html
// **********************************************************************

#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::Index;
using std::cout;
using std::endl;

/** Demo for initialization of dense matrices */

void initDense(Index const rows, Index const cols) {
  // Just allocate space for matrix, no initialisation; No guarantee for
  // zero-entries
  const MatrixXd A(rows, cols);
  // Zero matrix. Similar to matlab command zeros(rows, cols);
  const MatrixXd B = MatrixXd::Zero(rows, cols);
  // Ones matrix. Similar to matlab command ones(rows, cols);
  const MatrixXd C = MatrixXd::Ones(rows, cols);
  // Matrix with all entries same as value.
  const double value(3.14);
  const MatrixXd D = MatrixXd::Constant(rows, cols, value);
  // Random matrix, entries uniformly distributed in [0, 1]
  const MatrixXd E = MatrixXd::Random(rows, cols);
  // (Generalized) identity matrix, 1 on main diagonal
  const MatrixXd I = MatrixXd::Identity(rows, cols);

  cout << "size of A = (" << A.rows() << ',' << A.cols() << ')' << endl;
  cout << endl << "B = " << endl << B << endl;
  cout << endl << "C = " << endl << C << endl;
  cout << endl << "D = " << endl << D << endl;
  cout << endl << "E = " << endl << E << endl;
  cout << endl << "I = " << endl << I << endl;
}

int main() {
  cout << "initDense(3, 2)" << endl << endl;
  initDense(3, 2);
}