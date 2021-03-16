// **********************************************************************
// Eigen tutorial codes: http://eigen.tuxfamily.org/dox/GettingStarted.html
// **********************************************************************

#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/** Demo for use of matrix and vector types in Eigen */
/* SAM_LISTING_BEGIN_1 */
template <typename Scalar>
void eigenTypeDemo(unsigned int const dim) {
  // General dynamic (variable size) matrices
  using dynMat_t = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  // Dynamic (variable size) column vectors
  using dynColVec_t = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  // Dynamic (variable size) row vectors
  using dynRowVec_t = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;
  using index_t = typename dynMat_t::Index;
  using entry_t = typename dynMat_t::Scalar;

  // Declare vectors of size ´dim´; not yet initialized
  dynColVec_t colvec(dim);
  dynRowVec_t rowvec(dim);

  // Initialisation through component access
  for (index_t i = 0; i < colvec.size(); ++i) colvec(i) = (entry_t)i;

  for (index_t i = 0; i < rowvec.size(); ++i) rowvec(i) = (entry_t)1 / (i + 1);

  colvec[0] = (entry_t)3.14;
  rowvec[dim - 1] = (entry_t)2.718;

  // Form tensor product, a matrix, see Section 1.3.1
  dynMat_t vecprod = colvec * rowvec;
  const int nrows = vecprod.rows();
  const int ncols = vecprod.cols();

  cout << "colvec = " << endl
       << colvec << endl
       << endl
       << "rowvec = " << rowvec << endl
       << endl;

  cout << "size of vecprod = (" << nrows << ',' << ncols << ")" << endl;
}
/* SAM_LISTING_END_1 */

int main(int argc, char **argv) {
  cout << "eigenTypeDemo<float>(7)" << endl << endl;
  eigenTypeDemo<float>(7);
}