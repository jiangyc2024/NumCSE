// **********************************************************************
// Eigen tutorial codes: http://eigen.tuxfamily.org/dox/GettingStarted.html
// **********************************************************************

#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/** Demo for accessing sub-matrices */
/* SAM_LISTING_BEGIN_3 */
void accessColRow(int nrows = 6, int ncols = 7) {
  using index_t = typename Eigen::MatrixXd::Index;
  // Allocate dynamic \eigen matrix of type double
  MatrixXd m(nrows, ncols);
  // Initialization by direct component access
  for (index_t i = 0; i < m.rows(); i++)
    for (index_t j = 0; j < m.cols(); j++) m(i, j) = i + j;
  // Print matrix to standard output
  cout << "Matrix m = " << endl << m << endl << endl;
  // Print rows and columns
  for (index_t l = 0; l < m.rows(); l++)
    cout << "Row " << l << " = " << m.row(l) << endl;
  cout << endl;
  for (index_t l = 0; l < m.cols(); l++)
    cout << "Col " << l << " = " << endl << m.col(l) << endl << endl;
  // Access rows and columns as vectors
  RowVectorXd row1(m.row(1));
  VectorXd col1(m.col(1));

  cout << "Tensor product m.col(1) * m.row(1) = " << endl
       << col1 * row1 << endl
       << endl;

  m.col(2).setZero();
  m.row(2) = RowVectorXd::Constant(ncols, -1);
  m.row(4).tail(3) = RowVectorXd::Constant(3, 1.5);
  cout << "Modified matrix m = " << endl << m << endl;
}
/* SAM_LISTING_END_3 */

int main(int argc, char **argv) {
  cout << "accessColRow(6,7)" << endl << endl;
  accessColRow(6, 7);
}