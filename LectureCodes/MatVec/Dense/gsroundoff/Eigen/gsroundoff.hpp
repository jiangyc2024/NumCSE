#pragma once
#include <ios>
#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

#include "gramschmidt.hpp"

using namespace Eigen;
using namespace std;


//! \eigen{} Function demonstrating the effect of roundoff on the result
//! of Gram-Schmidt orthogonalization
//! A is 10x10 special matrix the so-called \href{https://en.wikipedia.org/wiki/Hilbert_matrix}{Hilbert matrix}:
//! $\MAc{i}{j} = (i+j-1)^{-1}$
/* SAM_LISTING_BEGIN_0 */
void gsroundoff(MatrixXd& A){
  // Gram-Schmidt orthogonalization of columns of A, see \cref{cpp:gramschmidt}
  MatrixXd Q = gramschmidt(A); 
  // \Magenta{Test orthonormality} of columns of Q, which should be an
  // \cor{orthogonal} matrix according to theory
  cout << setprecision(4) << fixed << "I = " 
       << endl << Q*Q.transpose() << endl; 
  // \eigen's \textbf{stable} internal Gram-Schmidt orthogonalization by
  // \cor{QR-decomposition}, see \cref{rem:QR} below
  HouseholderQR<MatrixXd> qr(A.rows(),A.cols());
  qr.compute(A); MatrixXd Q1 = qr.householderQ();
  // Test orthonormality
  cout << "I1 = " << endl << Q1*Q1.transpose() << endl;
  // Check orthonormality and span property \eqref{gsorth:span}
  MatrixXd R1 = qr.matrixQR().triangularView<Upper>();
  cout << scientific << "A-Q1*R1 = " << endl << A-Q1*R1 << endl;
}
/* SAM_LISTING_END_0 */
