/* **********************************************************************
 * Homework problem code for NumCSE course @ ETHZ
 * Lecturer: Prof. R. Hiptmair
 */
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
  // Efficient associative matrix multiplication
  MX = _A * (_B.transpose() * X);
  // END
  return MX;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
MatrixLowRank &MatrixLowRank::operator*=(const Eigen::MatrixXd &X) {
  // TO DO: (0-1.d)
  // START
  _B = X.transpose() * _B;
  _n = X.cols();
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
  unsigned int r = _r, s = X._r;
  // We will have X+this = Atilde*Btilde.transpose()
  MatrixXd Atilde(_m, r + s), Btilde(_n, r + s);
  Atilde << _A, X._A;
  Btilde << _B, X._B;

  // Calculate the QR decomposition of Atilde and Btilde.
  HouseholderQR<MatrixXd> qrA(Atilde);
  HouseholderQR<MatrixXd> qrB(Btilde);
  MatrixXd RA = qrA.matrixQR().triangularView<Upper>();
  MatrixXd RB = qrB.matrixQR().triangularView<Upper>();

  // We need the full SVD of RR = RA * RB.transpose()
  JacobiSVD<MatrixXd> svdRR(RA * RB.transpose(), ComputeFullU | ComputeFullV);

  // Rank of the sum X+this is the (numerical) rank of RA*RB.transpose().
  VectorXd Sigma = svdRR.singularValues();
  double sigma1 = Sigma[0];
  // Take care of the case sigma1=0, i.e. X+this = 0.
  // Follow the convention for storing the zero matrix
  if (sigma1 < atol) {
    _r = 1;
    _A = MatrixXd::Zero(_m, 1);
    _B = MatrixXd::Zero(_n, 1);
    return *this;
  }

  unsigned int rtilde = 1;
  // Numerical rank is at least 1.
  while (rtilde < Sigma.size() && Sigma[rtilde] > atol &&
         Sigma[rtilde] > rtol * sigma1) {
    // Numerical rank is at least rtilde+1.
    rtilde++;
  }

  _r = rtilde;
  _A = qrA.householderQ() * svdRR.matrixU().leftCols(rtilde) *
       Sigma.head(rtilde).asDiagonal();
  _B = qrB.householderQ() * svdRR.matrixV().leftCols(rtilde);
  // END
  return *this;
}
/* SAM_LISTING_END_5 */
