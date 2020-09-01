#include <iostream>
#include <algorithm>
#include <Eigen/Dense>
#include "matrixlowrank.hpp"

using namespace Eigen;

int main() {

  double tol = 1e-4;
  MatrixXd A(4,2),B(3,2),C(3,2),D(4,2),E(3,1),F(4,1),G(3,4);
  A << 1,2,3,4,5,6,7,8;
  B << 9,0,1,2,3,4;
  C << 5,6,7,8,9,0;
  D << 179,94,437,250,695,406,953,562;
  E << 9.+1e-9,1+1e-9,3+1e-9;
  F << -1,-3,-5,-7;  
  G << 0,0,0,0,4,8,12,16,8,16,24,32;
  MatrixLowRank L(B,A);
  MatrixLowRank M(A,B);
  MatrixLowRank N(E,F);
  
  /*
   *  run operator*()
   */
  MatrixXd MC = M*C;
  std::cout << "M*C = \n" << MC << std::endl;
  /*
   *  test operator*()
   */
  if (MC.cols() == D.cols() && MC.rows() == D.rows()) {
    if ( (MC-D).norm() < tol )
      std::cout << "Test for operator*() passed!\n\n";
    else
      std::cout
          << "Test for operator*() failed: wrong output.\n\n";
  } else
    std::cout
        << "Test for operator*() failed: wrong size.\n\n";

  /*
   *  run operator*=()
   */
  M *= C;
  MatrixXd Me = M*MatrixXd::Identity(M.cols(),M.cols()); // Eigen version of M.
  std::cout << "M after M*=C:\n" << Me << std::endl;
  /*
   *  test operator*=()
   */
  // This test can only be passed if 'operator*' is correct.
  if (Me.cols() == D.cols() && Me.rows() == D.rows()) {
    if ( (Me-D).norm() < tol && M.cols() == D.cols() && M.rows() == D.rows())
      std::cout << "Test for operator*=() passed!\n\n";
    else
      std::cout
          << "Test for operator*=() failed.\n\n";
  } else
    std::cout
        << "Test for operator*=() failed: wrong size.\n\n";

  /*
   *  run addTo()
   */
  N.addTo(L);
  MatrixXd Ne = N*MatrixXd::Identity(N.cols(),N.cols()); // Eigen version of N
  std::cout << "N after N.addTo(L):\n" << Ne << std::endl;
  JacobiSVD<MatrixXd> svd(G, ComputeFullU | ComputeFullV);
  svd.setThreshold(1e-6);
  /*
   *  test addTo()
   */
  // This test can only be passed if 'operator*' is correct.
  if (Ne.cols() == G.cols() && Ne.rows() == G.rows()) {
    if ( (Ne - G).norm() < tol && N.rank() == svd.rank() && N.cols() == G.cols() && N.rows() == G.rows())
      std::cout << "Test for addTo() passed!\n\n";
    else
      std::cout
          << "Test for addTo() failed.\n\n";
  } else
    std::cout
        << "Test for addTo() failed: wrong size.\n\n";  

  /*
   */   
  return 0;
}
