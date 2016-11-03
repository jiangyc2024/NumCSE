///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>

#include <Eigen/Dense>
#include <figure/figure.hpp>

#include "decomp.hpp"
#include "timer.h"

using namespace Eigen;

//! Eigen script for timing a smart and foolish way to carry out
//! multiplication with a scaling matrix
void housetime() {
/* SAM_LISTING_BEGIN_0 */
  int nruns = 3, minExp = 2, maxExp = 6;
  MatrixXd tms(maxExp-minExp+1,4);
  for(int i = 0; i <= maxExp-minExp; ++i){
    Timer t1, t2, t3;	// timer class
    int n = std::pow(2, minExp + i); int m = n*n;
    // Initialization of matrix A
    MatrixXd A(m,n); A.setZero();
    A.setIdentity(); A.block(n,0,m-n,n).setOnes();
    A += VectorXd::LinSpaced(m,1,m) * RowVectorXd::LinSpaced(n,1,n);
    for(int j = 0; j < nruns; ++j) {
      // plain QR-factorization in the constructor
      t1.start(); HouseholderQR<MatrixXd> qr(A); t1.stop();
      // full decomposition
      t2.start(); std::pair<MatrixXd,MatrixXd> QR2 = qr_decomp_full(A); t2.stop();
      // economic decomposition
      t3.start(); std::pair<MatrixXd,MatrixXd> QR3 = qr_decomp_eco(A); t3.stop();
    }
    tms(i,0)=n;
    tms(i,1)=t1.min(); tms(i,2)=t2.min(); tms(i,3)= t3.min(); 
  }
/* SAM_LISTING_END_0 */
  std::cout << std::scientific << std::setprecision(3) << tms << std::endl;

  // Plotting times using the mathGL and the Figure class.
  mgl::Figure fig;
  fig.setFontSize(4);
  fig.title("Tms for different ways to do qr decompositions");
  fig.setlog(true, true);
  fig.plot(tms.col(0),tms.col(1), " *r").label("HouseholderQR");
  fig.plot(tms.col(0),tms.col(2)," +b").label("qr decomposition full");
  fig.plot(tms.col(0),tms.col(3)," dm").label("qr decompositon eco");
  fig.xlabel("n");  fig.ylabel("time [s]");
  fig.legend(0.05,0.95);  fig.save("housetime_cpp");
}

