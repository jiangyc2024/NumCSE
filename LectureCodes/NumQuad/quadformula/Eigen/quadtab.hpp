///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::Matrix;

/* SAM_LISTING_BEGIN_0 */
class QuadTable {
public:
  using idx_t = MatrixXd::Index;
  const idx_t nmax = 10; // Maximal number of quadrature points 
private:
  const 
  static Matrix<double,nmax,nmax> QN; // nodes
  static Matrix<double,nmax,nmax> QW; // weights

  /*
http://stackoverflow.com/questions/31549398/c-eigen-initialize-static-matrix  static Eigen::Matrix4d foo = [] {
    Eigen::Matrix4d tmp;
    tmp << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
    return tmp;
}();
  */


/* SAM_LISTING_END_0 */
