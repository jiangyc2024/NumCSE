// **********************************************************************
// Eigen tutorial codes: http://eigen.tuxfamily.org/dox/GettingStarted.html
// **********************************************************************

#include <Eigen/Dense>
#include <iostream>

using Eigen::Matrix;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Dynamic;
using Eigen::ColMajor;
using Eigen::RowMajor;
using Eigen::Map;

using std::cout;
using std::endl;

/** Explores internal matrix storage format of Eigen */
/* SAM_LISTING_BEGIN_5 */
void storageOrder(int nrows = 6, int ncols = 7) {
  // Template parameter \texttt{ColMajor} selects column major data layout
  Matrix<double, Dynamic, Dynamic, ColMajor> mcm(nrows, ncols);
  // Template parameter \texttt{RowMajor} selects row major data layout
  Matrix<double, Dynamic, Dynamic, RowMajor> mrm(nrows, ncols);
  // Direct initialization; lazy option: use \texttt{int} as index type
  for (int l = 1, i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++, l++) {
      mcm(i, j) = mrm(i, j) = l;
    }
  }

  cout << "Matrix mrm = " << endl << mrm << endl << endl;
  cout << "mcm linear = ";
  for (int l = 0; l < mcm.size(); l++) {
    cout << mcm(l) << ' ';
  }
  cout << endl << endl;

  cout << "mrm linear = ";
  for (int l = 0; l < mrm.size(); l++) {
    cout << mrm(l) << ' ';
  }
  cout << endl << endl;

  // Retrieve pointer to raw matrix data
  double *mdat = mcm.data();
  cout << "mcm raw data layout in memory = ";
  for (int l = 0; l < mcm.size(); l++) {
    cout << mdat[l] << ' '; //NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  }
  cout << endl << endl;

  // Creates a map referencing the raw data of mcm
  // Reshaping possible only for matrices with non-prime dimensions
  Map<MatrixXd> mrs(mdat, nrows / 2, ncols * 2);
  cout << "mcm reshaped (mrs) = " << endl << mrs << endl << endl;
  // (Deep) copy data of mcm into matrix of different size (by implicit cast)
  MatrixXd cpy = mrs;
  mrs *= -1.5;
  cout << "mrs * -1.5 = " << endl << mrs << endl << endl;
  // Modifying mrs affects mcm, because they share the data space
  cout << "mcm modified = " << endl << mcm << endl << endl;
  // Matrix cpy is not affected, because of deep copy
  cout << "deep copy (cpy) = " << endl << cpy << endl;
}
/* SAM_LISTING_END_5 */

int main() {
  cout << "storageOrder(6, 7)" << endl;
  cout << "Demonstrating different storage orders for Eigen matrices" << endl
       << endl;
  storageOrder(6, 7);
}