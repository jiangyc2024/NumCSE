#pragma once

#include <iostream>

#include <cmath>

#include <Eigen/Dense>
//#include <figure/figure.hpp>
#include <mgl2/mgl.h>

#include "timer.h"

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! Eigen script for timing a smart and foolish way to carry out
//! multiplication with a scaling matrix
void scaletiming4() {
  int nruns = 3, minExp = 2, maxExp = 4;
  MatrixXd timings(maxExp-minExp+1,4);
  for(int i = 0; i <= maxExp-minExp; ++i){
    Timer tbad, tgood, topt;	// timer class
    int n = std::pow(2, minExp + i);
    VectorXd d = VectorXd::Random(n,1), x = VectorXd::Random(n,1), y(n);
    for(int j = 0; j < nruns; ++j) {
      MatrixXd D = d.asDiagonal(); // \label[Line]{scti:1}
      // matrix vector multiplication 
      tbad.start();  y = D*x; tbad.stop(); // \label[Line]{scti:2}
      // componentwise multiplication
      tgood.start(); y= d.cwiseProduct(x); tgood.stop(); // \label[Line]{scti:3}
      // matrix multiplication optimized by Eigen
      topt.start();  y = d.asDiagonal()*x; topt.stop(); // \label[Line]{scti:4}
    }
    timings(i,0)=n;
    timings(i,1)=tgood.min(); timings(i,2)=tbad.min(); timings(i,3)= topt.min(); 
  }
  std::cout << timings << std::endl;
  // Plotting times using the mathGL 
  mglGraph gr;
  //gr.LoadFont("heros","./Fonts/");
  gr.LoadFont("none");
  gr.LoadFont("heros");
  gr.SetFontSizePT(6);

  VectorXd x = timings.col(0),
           t1 = timings.col(1),
           t2 = timings.col(2),
           t3 = timings.col(3);

  mglData xd(x.data(), x.size()),
          t1d(t1.data(), t1.size()),
          t2d(t2.data(), t2.size()),
          t3d(t2.data(), t3.size());

  const double min = std::min( t1d.Minimal(), std::min( t2d.Minimal(), t3d.Minimal() ) ),
               max = std::max( t1d.Maximal(), std::max( t2d.Maximal(), t3d.Maximal() ) );

  gr.SetRanges( xd.Minimal(), xd.Maximal(), min, max );
  gr.SetFunc( "lg(x)", "lg(y)" );
  gr.Box();
  gr.Axis();
  gr.Title("Timings");
  gr.Plot(xd, t1d, " *r");
  gr.Plot(xd, t2d, " +b");
  gr.Plot(xd, t3d, " dm");
  gr.AddLegend("D.diagonal().cwiseProduct(x)", " *r");
  gr.AddLegend("D*x", " +b");
  gr.AddLegend("d.asDiagonal()*x", "dm");
  gr.Legend();
  gr.WriteEPS("sc4.eps");
              
  
  /*
  mgl::Figure fig;
  fig.setFontSize(4);
  fig.title("Timings for different ways to do scaling");
  fig.setlog(true, true);
  fig.plot(timings.col(0),timings.col(1), " *r").label("D.diagonal().cwiseProduct(x)");
  fig.plot(timings.col(0),timings.col(2)," +b").label("D*x");
  fig.plot(timings.col(0),timings.col(3)," dm").label("d.asDiagonal()*x");
  fig.xlabel("vector length n");  fig.ylabel("time [s]");
  fig.legend(0.05,0.95);  fig.save("sc0");
  */
}
/* SAM_LISTING_END_0 */
