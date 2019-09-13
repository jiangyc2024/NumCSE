///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <cmath>

#include <Eigen/Dense>
#include <figure/figure.hpp>

using namespace Eigen;

void linesec() {

  /* SAM_LISTING_BEGIN_0 */
  VectorXd phi = VectorXd::LinSpaced(50, M_PI / 200, M_PI / 2);
  MatrixXd res(phi.size(), 3);
  Matrix2d A;
  A(0, 0) = 1;
  A(1, 0) = 0;
  for (int i = 0; i < phi.size(); ++i) {
    A(0, 1) = std::cos(phi(i));
    A(1, 1) = std::sin(phi(i));
    // L2 condition number is the quotient of the maximal
    // and minimal singular value of A
    JacobiSVD<MatrixXd> svd(A);
    double C2 = svd.singularValues()(0) / // \Label[line]{cmc:1}
                svd.singularValues()(svd.singularValues().size() - 1);
    // L-infinity condition number
    double Cinf =
        A.inverse().cwiseAbs().rowwise().sum().maxCoeff() *
        A.cwiseAbs().rowwise().sum().maxCoeff(); // \Label[line]{cmc:2}
    res(i, 0) = phi(i);
    res(i, 1) = C2;
    res(i, 2) = Cinf;
  }
  /* SAM_LISTING_END_0 */
  // Plot
  mgl::Figure fig;
  fig.plot(res.col(0), res.col(1), "r").label("2-norm");
  fig.plot(res.col(0), res.col(2), "=b").label("max-norm");
  fig.xlabel("angle of n_1, n_2");
  fig.ylabel("condition numbers");
  fig.legend(1, 1);
  fig.save("linsec");
}
