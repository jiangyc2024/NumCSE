/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: December 2022
 */

#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <tuple>
#include <unsupported/Eigen/FFT>
#include <vector>

namespace SVDCirculant {
using namespace std::complex_literals;

/* SAM_LISTING_BEGIN_1 */
std::vector<unsigned int> permrevsort(Eigen::VectorXd &v) {
  // Retrieve length of vector
  const int n = v.size();
  // vector of integers 0, ... , n -1: stores permutation
  std::vector<unsigned int> perm(n);
  for (int j = 0; j < n; ++j) {
    perm[j] = j;
  }
  std::sort(perm.begin(), perm.end(),
            [&v](const unsigned int &i1, const unsigned int &i2) -> bool {
              return !(v[i1] < v[i2]);
            });
  const Eigen::VectorXd tmp{v};
  for (int j = 0; j < n; ++j) {
    v[j] = tmp[perm[j]];
  }
  return perm;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd, Eigen::MatrixXcd> svdCirculant(
    const Eigen::VectorXcd &u) {
  int n = u.size();  // Size of vectors and matrices
  Eigen::MatrixXcd U(n, n);
  Eigen::MatrixXcd V(n, n);
  // Apply \lref{lem:zirkdiag}: Compute entries of diagonal matrix via FFT
  Eigen::FFT<double> fft;
  Eigen::VectorXcd d_cplx = fft.fwd(u);
  // Diagonal of diagonal matrix in SVD (singular values)
  Eigen::VectorXd d = d_cplx.cwiseAbs();
  // Sort singular values
  auto perm = permrevsort(d);
  // Initialize Fourier matrix
  std::complex<double> omega =
      std::exp((2.0 * M_PI * 1.0i) / static_cast<double>(n));
  std::complex<double> fac = std::complex(1.0, 0.0);
  for (int k = 0; k < n; ++k) {
    std::complex<double> base =
        std::complex(1.0, 0.0) / std::sqrt(static_cast<double>(n));
    std::complex<double> d_unit = d_cplx[k] / std::abs(d_cplx[k]);
    for (int l = 0; l < n; ++l) {
      U(l, k) = base * d_unit;
      V(k, perm[l]) = std::conj(base);
      base *= fac;
    }
    fac *= omega;
  }
  return {U, d, V};
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd, Eigen::MatrixXcd> psvdCirculant(
    const Eigen::VectorXcd &u) {
  int n = u.size();          // Size of vectors and matrices
  Eigen::MatrixXcd U(n, n);  // U factor
  Eigen::MatrixXcd V(n, n);  // V factor
  // Apply \lref{lem:zirkdiag}: Compute entries of diagonal matrix via FFT
  Eigen::FFT<double> fft;
  Eigen::VectorXcd d_cplx = fft.fwd(u);
  // Diagonal of diagonal matrix in proto-SVD (singular values)
  Eigen::VectorXd d = d_cplx.cwiseAbs();
  // Initialize (scaled) Fourier matrix
  std::complex<double> omega =
      std::exp((2.0 * M_PI * 1.0i) / static_cast<double>(n));
  std::complex<double> fac = std::complex(1.0, 0.0);
  for (int k = 0; k < n; ++k) { // \Label[line]{psvdc:1}
    std::complex<double> base =
        std::complex(1.0, 0.0) / std::sqrt(static_cast<double>(n));
    // Diagonal entry of $\cob{\VD_\varphi}$
    std::complex<double> d_unit = d_cplx[k] / std::abs(d_cplx[k]);
    for (int l = 0; l < n; ++l) {
      U(l, k) = base * d_unit;
      V(k, l) = std::conj(base);
      base *= fac;
    }
    fac *= omega;
  } // \Label[line]{psvdc:2}
  return {U, d, V};
}
/* SAM_LISTING_END_3 */
/* SAM_LISTING_BEGIN_4 */
Eigen::MatrixXcd buildCircMat(const Eigen::VectorXcd &u) {
  int n = u.size();
  Eigen::MatrixXcd C(n, n);
  int j = 0;
  for (int k = 0; k < n; ++k, j = k) {
    for (int l = 0; l < n; ++l) {
      C(l, k) = u[j++];
      if (j >= n) {
        j = 0;
      }
    }
  }
  return C;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
bool testSvdCirculant(const Eigen::VectorXcd &u, double tol = 1.0E-10) {
  int n = u.size();
  auto [U, d, V] = psvdCirculant(u);
  // Test matrix factorization
  const Eigen::MatrixXcd C{buildCircMat(u)};
  if ((C - U * d.asDiagonal() * (V.conjugate().transpose())).norm() >
      tol * n * u.norm()) {
    std::cerr << "Factorization failed!" << std::endl;
    return false;
  }
  // Test whether matrices are unitary
  if ((Eigen::MatrixXcd::Identity(n, n) - U * U.conjugate().transpose())
          .norm() > tol) {
    std::cerr << "Matrix U is not unitary!" << std::endl;
    return false;
  }
  if ((Eigen::MatrixXcd::Identity(n, n) - V * V.conjugate().transpose())
          .norm() > tol) {
    std::cerr << "Matrix V is not unitary!" << std::endl;
    return false;
  }
  // Test whether singular values are positive
  for (int l = 0; l < n; ++l) {
    /*
      if (d[l] < d[l + 1]) {
      std::cerr << "Singular values not sorted" << std::endl;
      return false;
      }
      } */
    if (d[l] < 0.0) {
      std::cerr << "Negative singular values" << std::endl;
      return false;
    }
  }
  return true;
}
/* SAM_LISTING_END_5 */

}  // namespace SVDCirculant

int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "SVD of a circulant matrix" << std::endl;

  std::cout << "## Test of sorting function" << std::endl;
  Eigen::VectorXd v(5);
  v << 4.0, 3.0, 2.0, 5.0, 1.0;
  auto perm = SVDCirculant::permrevsort(v);
  std::cout << "Sorted vector v= " << v.transpose() << std::endl;
  std::cout << "Permutation: ";
  for (unsigned int &i : perm) {
    std::cout << i << " ";
  }
  std::cout << std::endl;

  std::cout << "## Test of building of circulant matrix" << std::endl;
  Eigen::VectorXcd u(5);
  u << 1.0, 2.0, 3.0, 4.0, 5.0;
  std::cout << "C(u) = " << std::endl
            << SVDCirculant::buildCircMat(u) << std::endl;

  std::cout << "## Test of SVD of circulant matrix" << std::endl;
  std::cout << "Test "
            << (SVDCirculant::testSvdCirculant(u) ? "passed" : "failed")
            << std::endl;
  return 0;
}
