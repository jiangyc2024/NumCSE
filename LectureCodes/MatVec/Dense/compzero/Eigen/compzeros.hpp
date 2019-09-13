///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <cmath>
#include <iostream>

#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace std;
using namespace Eigen;

#include "zerosquadpol.hpp"
#include "zerosquadpolstab.hpp"

/* SAM_LISTING_BEGIN_0 */
//! Eigen Function for testing the computation of the zeros of a parabola
void compzeros() {
  int n = 100;
  MatrixXd res(n, 4);
  VectorXd gamma = VectorXd::LinSpaced(n, 2, 992);
  for (int i = 0; i < n; ++i) {
    double alpha = -(gamma(i) + 1. / gamma(i));
    double beta = 1.;
    Vector2d z1 = zerosquadpol(alpha, beta);
    Vector2d z2 = zerosquadpolstab(alpha, beta);
    double ztrue = 1. / gamma(i), z2true = gamma(i);
    res(i, 0) = gamma(i);
    res(i, 1) = std::abs((z1(0) - ztrue) / ztrue);
    res(i, 2) = std::abs((z2(0) - ztrue) / ztrue);
    res(i, 3) = std::abs((z1(1) - z2true) / z2true);
  }
  /* SAM_LISTING_END_0 */
  // Graphical output of relative error of roots computed by unstable
  // implementation
  mgl::Figure fig1;
  fig1.setFontSize(3);
  fig1.title("Roots of a parabola computed in an unstable manner");
  fig1.plot(res.col(0), res.col(1), " +r").label("small root");
  fig1.plot(res.col(0), res.col(3), " *b").label("large root");
  fig1.xlabel("\\gamma");
  fig1.ylabel("relative errors in \\xi_1, \\xi_2");
  fig1.legend(0.05, 0.95);
  fig1.save("zqperrinstab");
  // Graphical output of relative errors (comparison), small roots
  mgl::Figure fig2;
  fig2.title("Roundoff in the computation of zeros of a parabola");
  fig2.plot(res.col(0), res.col(1), " +r").label("unstable");
  fig2.plot(res.col(0), res.col(3), " *m").label("stable");
  fig2.xlabel("\\gamma");
  fig2.ylabel("relative errors in \\xi_1");
  fig2.legend(0.05, 0.95);
  fig2.save("zqperr");
}
