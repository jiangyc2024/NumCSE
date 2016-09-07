#pragma once

#include <iostream>

#include <cmath>

#include <Eigen/Dense>
#include <figure/figure.hpp>

#include "timer.h"

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! Eigen script for timing a smart and foolish way to carry out
//! multiplication with a scaling matrix
void scaletiming(){
  int nruns = 3, minExp = 2, maxExp = 14;
  MatrixXd timings(maxExp-minExp+1,3);
  for(int i = 0; i <= maxExp-minExp; ++i){
	Timer tbad, tgood;	// timer class
	int n = std::pow(2, minExp + i);
    VectorXd d = VectorXd::Random(n,1), x = VectorXd::Random(n,1), y(n);
    for(int j = 0; j < nruns; ++j){
		MatrixXd D = d.asDiagonal(); // d.asDiagonal()*x would be optimized \label{scti:1}
		tbad.start(); y = D*x; 	tbad.stop();	// matrix vector multiplication \label{scti:2}
		tgood.start();y= d.cwiseProduct(x);	tgood.stop(); // comp. wise mult. \label{scti:3}
	}
	timings(i,0)=n; timings(i,1)=tgood.min(); timings(i,2)=tbad.min();
  }
  std::cout << timings << std::endl;
  //Plotting
  mgl::Figure fig;
  fig.setFontSize(4);
  fig.title("Timings for different ways to do scaling");
  fig.setlog(true, true);
  fig.plot(timings.col(0),timings.col(1), " *r").label("D.diagonal().cwiseProduct(x)");
  fig.plot(timings.col(0),timings.col(2)," +b").label("D*x");
  fig.xlabel("vector length n");  fig.ylabel("time [s]");
  fig.legend(0.05,0.95);  fig.save("scaletiming");
}
/* SAM_LISTING_END_0 */
