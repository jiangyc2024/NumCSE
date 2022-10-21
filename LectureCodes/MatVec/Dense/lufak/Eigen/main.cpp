///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <iostream>
#include <utility>

#include "lufak.hpp"

/* SAM_LISTING_BEGIN_1 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
lufak_eigen(const Eigen::MatrixXd &A) {
  // Compute LU decomposition 
  auto ludec = A.lu();
  // The LU-factors are computed by in-situ LU-decomposition,
  // see \cref{rem:insitulu}, and are stored in a dense matrix of
  // the same size as A 
  Eigen::MatrixXd L { ludec.matrixLU().triangularView<Eigen::UnitLower>() };
  Eigen::MatrixXd U { ludec.matrixLU().triangularView<Eigen::Upper>() };
  // \eigen employs partial pivoting, see \cref{alg:GEp}, which can be viewed
  // as a prior permutation of the rows of A. We apply the inverse of this
  // permutation to the L-factor in order to achieve \cob{$\VA=\VL\VU$}. 
  L.applyOnTheLeft(ludec.permutationP().inverse());
  // Return LU-factors as members of a 2-tuple. 
  return { L , U }; 
}
/* SAM_LISTING_END_1 */

int main() {
  int n = 3;
  Eigen::MatrixXd A(n, n);
  A << 1, 0, 2, -1, 4, 1, -2, 1, 2;
  {
    std::cout << "(1) Hand-made LU factorization" << std::endl;
    auto [L,U] = lufak::lufak(A);
    std::cout << "L=\n" << L << std::endl << "U=\n" << U
	     << "|A-L*U| = " << (A-L*U).norm() << std::endl;
  }

  {
    std::cout << "(2) (Permuted) LU factors computed by Eigen" << std::endl;
    auto [L,U] = lufak_eigen(A);
    std::cout << "L=\n" << L << std::endl << "U=\n" << U
	     << "|A-L*U| = " << (A-L*U).norm() <<std::endl;
  }
  return 0;
}
