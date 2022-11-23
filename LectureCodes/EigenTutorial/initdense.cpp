// **********************************************************************
// Eigen tutorial codes: http://eigen.tuxfamily.org/dox/GettingStarted.html
// **********************************************************************

#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/** Demo for initialization of dense matrices */

void initDense(size_t const rows, size_t const cols) {
  // Just allocate space for matrix, no initialisation; No guarantee for
  // zero-entries
  MatrixXd A(rows, cols);
  // Zero matrix. Similar to matlab command zeros(rows, cols);
  MatrixXd B = MatrixXd::Zero(rows, cols);
  // Ones matrix. Similar to matlab command ones(rows, cols);
  MatrixXd C = MatrixXd::Ones(rows, cols);
  // Matrix with all entries same as value.
  double value(3.14);
  MatrixXd D = MatrixXd::Constant(rows, cols, value);
  // Random matrix, entries uniformly distributed in [0, 1]
  MatrixXd E = MatrixXd::Random(rows, cols);
  // (Generalized) identity matrix, 1 on main diagonal
  MatrixXd I = MatrixXd::Identity(rows, cols);

  cout << "size of A = (" << A.rows() << ',' << A.cols() << ')' << endl;
  cout << endl << "B = " << endl << B << endl;
  cout << endl << "C = " << endl << C << endl;
  cout << endl << "D = " << endl << D << endl;
  cout << endl << "E = " << endl << E << endl;
  cout << endl << "I = " << endl << I << endl;
}

int main(int argc, char **argv) {
  cout << "initDense(3, 2)" << endl << endl;
  initDense(3, 2);
}