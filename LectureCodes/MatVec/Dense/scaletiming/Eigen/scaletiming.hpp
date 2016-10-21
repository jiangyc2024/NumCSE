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

#include "timer.h"

using namespace Eigen;

//! Eigen script for timing a smart and foolish way to carry out
//! multiplication with a scaling matrix
void scaletiming() {
/* SAM_LISTING_BEGIN_0 */
  int nruns = 3, minExp = 2, maxExp = 14;
  MatrixXd tms(maxExp-minExp+1,4);
  for(int i = 0; i <= maxExp-minExp; ++i){
    Timer tbad, tgood, topt;	// timer class
    int n = std::pow(2, minExp + i);
    VectorXd d = VectorXd::Random(n,1), x = VectorXd::Random(n,1), y(n);
    for(int j = 0; j < nruns; ++j) {
      MatrixXd D = d.asDiagonal(); // \Label[Line]{scti:1}
      // matrix vector multiplication 
      tbad.start();  y = D*x; tbad.stop(); // \Label[Line]{scti:2}
      // componentwise multiplication
      tgood.start(); y= d.cwiseProduct(x); tgood.stop(); // \Label[Line]{scti:3}
      // matrix multiplication optimized by Eigen
      topt.start();  y = d.asDiagonal()*x; topt.stop(); // \Label[Line]{scti:4}
    }
    tms(i,0)=n;
    tms(i,1)=tgood.min(); tms(i,2)=tbad.min(); tms(i,3)= topt.min(); 
  }
/* SAM_LISTING_END_0 */
  std::cout << std::scientific << std::setprecision(3) << tms << std::endl;

  // Plotting times using the mathGL and the Figure class.
  mgl::Figure fig;
  fig.setFontSize(4);
  fig.title("Tms for different ways to do scaling");
  fig.setlog(true, true);
  fig.plot(tms.col(0),tms.col(1), " *r").label("D.diagonal().cwiseProduct(x)");
  fig.plot(tms.col(0),tms.col(2)," +b").label("D*x");
  fig.plot(tms.col(0),tms.col(3)," dm").label("d.asDiagonal()*x");
  fig.xlabel("vector length n");  fig.ylabel("time [s]");
  fig.legend(0.05,0.95);  fig.save("scaletiming");
}

