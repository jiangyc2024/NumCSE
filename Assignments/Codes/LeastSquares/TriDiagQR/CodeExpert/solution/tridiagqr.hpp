#include <Eigen/Dense>
#include <iostream>

/* SAM_LISTING_BEGIN_0 */
struct TriDiagonalMatrix {
  Eigen::Index n;     // Matrix size $n\times n$
  Eigen::VectorXd d;  // $n$-vector of diagonal entries
  Eigen::VectorXd l;  // $n-1$-vector, entries of first lower diagonal
  Eigen::VectorXd u;  // $n-1$-vector, entries of first upper diagonal

  TriDiagonalMatrix(const Eigen::VectorXd &d, const Eigen::VectorXd &l,
                    const Eigen::VectorXd &u)
      : n(d.size()), d(d), l(l), u(u) {
    assert((n - 1) == l.size() && (n - 1) == u.size());
  }
};
/* SAM_LISTING_END_0 */

/* @brief Returns +1 for x>=0 and -1 for x<0.
 */
inline double sign(double x) {
  if (x >= 0)
    return 1.0;
  else
    return -1.0;
}

// An alternative:
// https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

/* @brief Calculates the rotation matrix corresponding to rho.
 */
/* SAM_LISTING_BEGIN_9 */
inline Eigen::Matrix2d Givens(double rho) {
  double gamma, sigma;
  if (std::abs(rho) < 1.0) {
    sigma = 0.5 * std::sqrt(2.0) * rho;
    gamma = std::sqrt(1.0 - sigma * sigma);
  } else if (std::abs(rho) > 1.0) {
    gamma = 0.5 * std::sqrt(2.0) / rho;
    sigma = std::sqrt(1 - gamma * gamma);
  } else {
    gamma = 0.0;
    sigma = sign(rho);
  }
  Eigen::Matrix2d G;
  G << gamma, sigma, -sigma, gamma;
  return G;
}
/* SAM_LISTING_END_9 */

/* @brief Calculates the parameters of a Givens rotation for a.
 */
/* SAM_LISTING_BEGIN_1 */
std::tuple<double, double, double> compGivensRotation(Eigen::Vector2d a) {
  double rho, gamma, sigma;
  // TO DO (4-4.b): Compute the parameters of a rotation that when
  // applied to $a$ results in the vector $(\|a\|, 0)^T$ or $(-\|a\|, 0)^T$.
  // Hint: You can use the function sign() above when determining rho
  // from gamma and sigma.
  // START
  if (a(1) != 0.0) {
    double t;
    // Separate cases to avoid cancellation/overflow
    if (std::abs(a(1)) > std::abs(a(0))) {
      t = -a(0) / a(1);
      sigma = 1.0 / std::sqrt(1.0 + t * t);
      gamma = sigma * t;
    } else {             // Then abs(a(0)) >= abs(a(1)) > 0
      t = -a(1) / a(0);  // a(0) != 0
      gamma = 1.0 / std::sqrt(1.0 + t * t);
      sigma = gamma * t;
    }
  } else {
    // No rotation.
    gamma = 1.0;
    sigma = 0.0;
  }
  // Determine rho
  if (gamma == 0.0)
    rho = sign(sigma);
  else if (std::abs(gamma) > std::abs(sigma))
    rho = std::sqrt(2.0) * sign(gamma) * sigma;
  else
    rho = 0.5 * std::sqrt(2.0) * sign(sigma) / gamma;
  // END
  return std::make_tuple(rho, gamma, sigma);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
class TriDiagonalQR {
 public:
  explicit TriDiagonalQR(const TriDiagonalMatrix &A);

  template <typename VecType>
  Eigen::VectorXd applyQT(const VecType &x) const;
  template <typename VecType>
  Eigen::VectorXd solve(const VecType &b) const;
  // For debugging purposes: extract factors as dense matrices
  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> getQRFactors(void) const;

 private:
  Eigen::Index n;       // size of the square matrix
  Eigen::MatrixXd B;    // Three non-zero upper diagonals of R
  Eigen::VectorXd rho;  // Encoded Givens rotations
};
/* SAM_LISTING_END_2 */

/* @brief Computes the QR decomposition of the tridiagonal matrix A.
 */
/* SAM_LISTING_BEGIN_3 */
TriDiagonalQR::TriDiagonalQR(const TriDiagonalMatrix &A) : n(A.n) {
  rho = Eigen::VectorXd::Zero(n - 1);
  B = Eigen::MatrixXd::Zero(n, 3);
  B.col(0) = A.d;
  B.col(1).tail(n - 1) = A.l;
  B.col(2).head(n - 1) = A.u;
  // Currently we have for k=0,...,n-1:
  // B(k,0) = A(k,k), B(k,1) = A(k,k-1), B(k,2) = A(k,k+1).

  // A bit clumsy, but less error prone than juggling indices
  Eigen::MatrixXd Ablock(2, 3);
  double gamma, sigma;
  for (int k = 0; k < (n - 1); k++) {
    // The block of the current A that will be affected by the next rotation.
    Ablock << B(k, 0), B(k, 2), 0, B(k + 1, 1), B(k + 1, 0), B(k + 1, 2);
    // TO DO (4-4.c): Find the k-th Givens rotation and update rho and B.
    // START
    std::tuple<double, double, double> params =
        compGivensRotation(Ablock.col(0));
    rho[k] = std::get<0>(params);  // Store rotation parameter
    gamma = std::get<1>(params);
    sigma = std::get<2>(params);

    // Update B
    Eigen::Matrix2d G;
    G << gamma, -sigma, sigma, gamma;
    Ablock = G * Ablock;  // Rotate to eliminate Ablock(1,0)
    B.row(k) = Ablock.row(0);
    B(k + 1, 0) = Ablock(1, 1);
    B(k + 1, 2) = Ablock(1, 2);
    // END
  }
}
/* SAM_LISTING_END_3 */

/* @brief Computes Q.transpose()*x.
 */
/* SAM_LISTING_BEGIN_4 */
template <typename VecType>
Eigen::VectorXd TriDiagonalQR::applyQT(const VecType &x) const {
  assert(n == x.size() &&
         "Input must have size equal to the number of rows in Q.");
  Eigen::VectorXd y(n);
  // TO DO (4-4.d): Fill in the entries of y.
  // START
  for (int k = 0; k < n; k++) {
    y[k] = x[k];
  }
  Eigen::Matrix2d G;
  for (int k = 0; k < (n - 1); k++) {
    G = Givens(rho[k]).transpose();
    y.segment(k, 2) = G * y.segment(k, 2);
  }
  // END
  return y;
}
/* SAM_LISTING_END_4 */

/* @brief Solves A*x=b.
 */
/* SAM_LISTING_BEGIN_5 */
template <typename VecType>
Eigen::VectorXd TriDiagonalQR::solve(const VecType &b) const {
  assert(n == b.size() &&
         "Input must have size equal to the number of rows in Q.");
  Eigen::VectorXd x(n);
  // TO DO (4-4.e): Solve the linear system A*x = b if A is invertible,
  // otherwise throw an std::runtime_error.
  // START
  if ((B.col(0).array().abs() < 1E-10).any()) {
    throw std::runtime_error("A is not invertible!");
  }
  // Multiply by Q^T on both sides of Q*R*x=b:
  // R*x = Q^T*b = c
  Eigen::VectorXd c = applyQT(b);

  // Back-substitution
  int k = n - 1;
  x(k) = c(k) / B(k, 0);
  k--;
  x(k) = (c(k) - B(k, 1) * x(k + 1)) / B(k, 0);
  k--;
  while (k >= 0) {
    x(k) = (c(k) - B(k, 1) * x(k + 1) - B(k, 2) * x(k + 2)) / B(k, 0);
    k--;
  }
  // END
  return x;
}
/* SAM_LISTING_END_5 */

/* @brief Computes Q and R as dense matrices.
 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> TriDiagonalQR::getQRFactors()
    const {
  // This function is only for demonstration purposes,
  // and is not in any way optimal.
  Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(n, n);
  R.diagonal() = B.col(0);
  R.diagonal(1) = B.col(1).head(n - 1);
  R.diagonal(2) = B.col(2).head(n - 2);
  Eigen::MatrixXd G(n, n);
  for (int k = 0; k < (n - 1); k++) {
    G.setIdentity();
    G.block(k, k, 2, 2) = Givens(rho[k]);
    Q = Q * G;
  }

  return std::make_pair(Q, R);
}

/* SAM_LISTING_BEGIN_6 */
template <typename MatrixType>
unsigned int invit(const MatrixType &A, Eigen::VectorXd &x, double TOL = 1E-6,
                   unsigned int maxit = 100) {
  x.normalize();
  Eigen::VectorXd x_old(x.size());
  int it = 0;
  do {
    x_old = x;
    x = A.lu().solve(x_old);
    x.normalize();
  } while (((x - x_old).norm() > TOL) && (it++ < maxit));
  return it;
}
/* SAM_LISTING_END_6 */

/* SAM_LISTING_BEGIN_7 */
template <>
unsigned int invit<TriDiagonalMatrix>(const TriDiagonalMatrix &A,
                                      Eigen::VectorXd &x, double TOL,
                                      unsigned int maxit) {
  int it = 0;
  // TO DO (4-4.g): Specialize (efficiently) the template function invit()
  // for the TriDiagonalMatrix class.
  // START
  x.normalize();
  Eigen::VectorXd x_old(x.size());
  TriDiagonalQR Aqr(A);
  do {
    x_old = x;
    x = Aqr.solve(x_old);
    x.normalize();
  } while (((x - x_old).norm() > TOL) && (it++ < maxit));
  // END
  return it;
}
/* SAM_LISTING_END_7 */
