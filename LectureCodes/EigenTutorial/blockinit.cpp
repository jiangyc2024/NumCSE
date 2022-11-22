// **********************************************************************
// Eigen tutorial codes: http://eigen.tuxfamily.org/dox/GettingStarted.html
// **********************************************************************

#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;

using std::cout;
using std::endl;

/** Demo for block initialization of matrices */
/* SAM_LISTING_BEGIN_2 */
void blockInit(int size = 6) {
  // Make size an even number
  size = 2 * (size / 2);

  MatrixXd mat(size, size);
  // Initialize matrix with the comma initializer ´<<´, from top-left to
  // bottom-right
  mat << MatrixXd::Zero(size / 2, size / 2),
      MatrixXd::Identity(size / 2, size / 2),
      MatrixXd::Identity(size / 2, size / 2),
      MatrixXd::Zero(size / 2, size / 2);
  cout << "#1 " << endl << mat << endl << endl;

  // Set all matrix entries to their column-major order
  MatrixXd botRows(size / 2, size);
  for (int l = 0; l < botRows.size(); ++l) {
    botRows(l) = l;
  }

  // Blocks can have different shapes
  mat << MatrixXd::Zero(size / 2, size / 2),
      MatrixXd::Identity(size / 2, size / 2), botRows;
  cout << "#2 " << endl << mat << endl << endl;

  // Stack two matrices
  mat << MatrixXd::Constant(size, size - 4, 1.5),
      MatrixXd::Constant(size, 4, 3.5);
  cout << "#3 " << endl << mat << endl << endl;

  mat << MatrixXd::Constant(size - 2, size - 4, 1.5),  // top row, first block
      MatrixXd::Constant(size - 2, 3, 3.5),            // top row, second block
      MatrixXd::Constant(size - 2, 1, 7.5),            //  top row, third block
      MatrixXd::Constant(2, size - 2, 2.5),            // bottom row, left block
      MatrixXd::Constant(2, 2, 4.5);  // bottom row, right block

  cout << "#4 " << endl << mat << endl << endl;
}
/* SAM_LISTING_END_2 */

int main() {
  cout << "blockInit(6)" << endl << endl;
  blockInit(6);
}