#include <Eigen/Dense>
#include <algorithm>
#include <iostream>

/* SAM_LISTING_BEGIN_0 */
class MatrixLowRank {
 public:
  MatrixLowRank(unsigned int m, unsigned int n, unsigned int r);
  MatrixLowRank(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B);
  MatrixLowRank() = default;

  Eigen::Index rows() const { return _m; };
  Eigen::Index cols() const { return _n; };
  Eigen::Index rank() const { return _r; };

  Eigen::MatrixXd operator*(const Eigen::MatrixXd &) const;
  MatrixLowRank &operator*=(const Eigen::MatrixXd &);
  MatrixLowRank &addTo(const MatrixLowRank &, double rtol = 1E-6,
                       double atol = 1E-8);

 private:
  unsigned int _m;     // number of rows
  unsigned int _n;     // number of columns
  unsigned int _r;     // maximal rank, =1 for zero matrix
  Eigen::MatrixXd _A;  // factor matrix A
  Eigen::MatrixXd _B;  // factor matrix B
};

/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
MatrixLowRank::MatrixLowRank(unsigned int m, unsigned int n, unsigned int r)
    : _m(m), _n(n), _r(r) {
  // Some arbitrary choice: not really meaningful
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
  Eigen::MatrixXd MX =
      Eigen::MatrixXd::Zero(_m, _n);  // will contain the result.
  // TODO: (3-11.b) Implement the multiplication of the $m \times n$-matrix
  // $M = AB^T$ stored in this class with the given eigen matrix X.
  // Note: The implementation should be efficient.

  // START
  // Efficient associative matrix multiplication
  MX = _A * (_B.transpose() * X);
  // END

  return MX;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
MatrixLowRank &MatrixLowRank::operator*=(const Eigen::MatrixXd &X) {
  // TODO: (3-11.d) Implement the *= operator for the in-situ multiplication of
  // the $m \times n$-matrix $M = AB^T$ stored in this class with the given
  // eigen matrix X.

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
  // TODO: (3-11.e) Implement the function s.t. the $m \times n$-matrix stored
  // in this class is replaced by trunc_tol(M + X).
  // Note: X is not an eigen matrix but a matrix of type MatrixLowRank.
  // Note: We have this = _A*_B.transpose()
  // Note: We have X = X._A*X._B.transpose()
  // Note: This function should be efficient.

  // START
  const unsigned int r = _r, s = X._r;
  // We will have X+this = Atilde*Btilde.transpose()
  Eigen::MatrixXd Atilde(_m, r + s), Btilde(_n, r + s);
  Atilde << _A, X._A;
  Btilde << _B, X._B;

  // Calculate the QR decomposition of Atilde and Btilde.
  Eigen::HouseholderQR<Eigen::MatrixXd> qrA(Atilde);
  Eigen::HouseholderQR<Eigen::MatrixXd> qrB(Btilde);
  Eigen::MatrixXd RA = qrA.matrixQR().triangularView<Eigen::Upper>();
  Eigen::MatrixXd RB = qrB.matrixQR().triangularView<Eigen::Upper>();

  // We need the full SVD of RR = RA * RB.transpose()
  Eigen::JacobiSVD<Eigen::MatrixXd> svdRR(
      RA * RB.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);

  // Rank of the sum X+this is the (numerical) rank of RA*RB.transpose().
  Eigen::VectorXd Sigma = svdRR.singularValues();
  const double sigma1 = Sigma[0];
  // Take care of the case sigma1=0, i.e. X+this = 0.
  // Follow the convention for storing the zero matrix
  if (sigma1 < atol) {
    _r = 1;
    _A = Eigen::MatrixXd::Zero(_m, 1);
    _B = Eigen::MatrixXd::Zero(_n, 1);

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
