///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <cassert>

#include <Eigen/Dense>

#include "bandwidth.hpp"
#include "substenv.hpp"

using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! envelope aware recursive LU-factorization
//! of structurally symmetric matrix
void luenv(const MatrixXd &A, MatrixXd &L, MatrixXd &U) {
  int n = A.cols();
  assert(n == A.rows() && "A must be square");
  if (n == 1) {
    L.setIdentity();
    U = A;
  } else {
    VectorXi mr = rowbandwidth(A);  // = colbandwidth thanks to symmetry
    double gamma;
    MatrixXd L1(n - 1, n - 1), U1(n - 1, n - 1);
    luenv(A.topLeftCorner(n - 1, n - 1), L1, U1);
    VectorXd u = substenv(L1, A.col(n - 1).head(n - 1), mr);
    VectorXd l =
        substenv(U1.transpose(), A.row(n - 1).head(n - 1).transpose(), mr);
    if (mr(n - 1) > 0)
      gamma = A(n - 1, n - 1) - l.tail(mr(n - 1)).dot(u.tail(mr(n - 1)));
    else
      gamma = A(n - 1, n - 1);
    L.topLeftCorner(n - 1, n - 1) = L1;
    L.col(n - 1).setZero();
    L.row(n - 1).head(n - 1) = l.transpose();
    L(n - 1, n - 1) = 1;
    U.topLeftCorner(n - 1, n - 1) = U1;
    U.col(n - 1).head(n - 1) = u;
    U.row(n - 1).setZero();
    U(n - 1, n - 1) = gamma;
  }
}
/* SAM_LISTING_END_0 */
