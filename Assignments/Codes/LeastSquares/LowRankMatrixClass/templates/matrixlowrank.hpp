#include <Eigen/Dense>
#include <algorithm>
#include <iostream>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
class MatrixLowRank {
public:
  MatrixLowRank(unsigned int m, unsigned int n, unsigned int r);
  MatrixLowRank(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B);

  Eigen::Index rows() const { return _m; };
  Eigen::Index cols() const { return _n; };
  Eigen::Index rank() const { return _r; };

  Eigen::MatrixXd operator*(const Eigen::MatrixXd &)const;
  MatrixLowRank &operator*=(const Eigen::MatrixXd &);
  MatrixLowRank &addTo(const MatrixLowRank &, double rtol = 1E-6,
                       double atol = 1E-8);

private:
  unsigned int _m;    // no. of rows
  unsigned int _n;    // no. of columns
  unsigned int _r;    // maximal rank, =1 for zero matrix
  Eigen::MatrixXd _A; // factor matrix A
  Eigen::MatrixXd _B; // factor matrix B
};

/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
MatrixLowRank::MatrixLowRank(unsigned int m, unsigned int n, unsigned int r)
    : _m(m), _n(n), _r(r) {
  // Some abritrary choice: not really meaningful
  _A = Eigen::MatrixXd::Identity(m, r);
  _B = Eigen::MatrixXd::Identity(n, r);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
MatrixLowRank::MatrixLowRank(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B)
    : _m(A.rows()), _n(B.rows()), _r(A.cols()), _A(A), _B(B) {
  assert(_r == B.cols() &&
         "No. of columns in A must be equal to no. of columns in B!");
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd MatrixLowRank::operator*(const Eigen::MatrixXd &X) const {
  MatrixXd MX; // will contain the result.
  // TO DO: (0-1.b)
  // START
  
  // END
  return MX;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
MatrixLowRank &MatrixLowRank::operator*=(const Eigen::MatrixXd &X) {
  // TO DO: (0-1.d)
  // START

  // END
  return *this;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
MatrixLowRank &MatrixLowRank::addTo(const MatrixLowRank &X, double rtol,
                                    double atol) {
  assert(_m == X._m && _n == X._n && "Matrices must have the same shape!");
  assert((rank() + X.rank() <= std::min(_n, _m)) && "Sum of ranks too large!");
  assert(rtol > 0 && atol > 0 && "Tolerances must be positive!");
  // TO DO: (0-1.e)
  // We have this = _A*_B.transpose() and
  // X = X._A*X._B.transpose()
  // START

  // END
  return *this;
}
/* SAM_LISTING_END_5 */
