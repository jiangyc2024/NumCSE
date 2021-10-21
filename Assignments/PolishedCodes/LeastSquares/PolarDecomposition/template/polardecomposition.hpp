#ifndef POLARDECOMPOSITIONHPP
#define POLARDECOMPOSITIONHPP

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <cassert>
#include <cmath>
#include <iostream>

class PolarDecomposition {
public:
  explicit PolarDecomposition(const Eigen::MatrixXd &X) { initialize(X); }
  PolarDecomposition(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B);
  PolarDecomposition(const PolarDecomposition &) = default;
  ~PolarDecomposition() = default;

  // Left multiplication of M with the Q-factor of the polar decomposition
  void applyQ(Eigen::MatrixXd &Y) const { Y.applyOnTheLeft(Q_); }
  // Left multiplication of M with the M-factor of the polar decomposition
  void applyM(Eigen::MatrixXd &Y) const { Y.applyOnTheLeft(M_); }
  
  int Qcols() { return Q_.cols(); }
  int Mcols() { return M_.cols(); }

private:
  void initialize(const Eigen::MatrixXd &X);
  Eigen::MatrixXd Q_; // factor Q
  Eigen::MatrixXd M_; // factor M
};

/* SAM_LISTING_BEGIN_1 */
void PolarDecomposition::initialize(const Eigen::MatrixXd &X) {
  assert(X.rows() >= X.cols());
  // TO DO: Implement a method to initialize the data members _Q and _M 
  // corresponding to Q and M in Theorem 0.3.1, where X = QM
  // START
  
  // END
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_7 */
PolarDecomposition::PolarDecomposition(const Eigen::MatrixXd &A,
                                       const Eigen::MatrixXd &B) {
  const long m = A.rows();  // No. of rows of $\cob{\VX}$
  const long n = B.rows();  // No. of columns of $\cob{\VX}$
  const long k = A.cols();  // Maximal rank of $\cob{\VX}$
  // We assume $\cob{k \leq n \leq m}$
  assert(m >= n);
  assert(k < n);
  assert(B.cols() == k);
  // TO DO: Implement a method to initialize the data members _Q and _M
  // for X := AB^T = QM, with optimal complexity
  // START
  
  // END
}
/* SAM_LISTING_END_7 */


#endif
