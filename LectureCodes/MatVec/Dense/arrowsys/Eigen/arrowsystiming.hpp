///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include "arrowsys_fast.hpp"
#include "arrowsys_slow.hpp"
#include "arrowsys_sparse.hpp"
#include "timer.h"

namespace arrowsys {


using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseLU;
using Eigen::SparseMatrix;
using Eigen::BiCGSTAB;

inline
/* SAM_LISTING_BEGIN_0 */
MatrixXd arrowsystiming() {
  std::vector<int> n = {8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
  const int nruns = 3;
  MatrixXd times(n.size(), 6);
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(n.size()); ++i) {
    // timer class
    Timer t1;
    Timer t2;
    Timer t3;
    Timer t4; 
    const double alpha = 2;
    const VectorXd b = VectorXd::Ones(n[i], 1);
    const VectorXd c = VectorXd::LinSpaced(n[i], 1, n[i]);
    const VectorXd d = -b;
    const VectorXd y = VectorXd::Constant(n[i] + 1, -1)
                     .binaryExpr(VectorXd::LinSpaced(n[i] + 1, 1, n[i] + 1),
                                 [](double x, double y) { return pow(x, y); })
                     .array();
    VectorXd x1(n[i] + 1);
    VectorXd x2(n[i] + 1);
    VectorXd x3(n[i] + 1);
    VectorXd x4(n[i] + 1);
    for (int j = 0; j < nruns; ++j) {
      t1.start();
      x1 = arrowsys_slow(d, c, b, alpha, y);
      t2.stop();
      t2.start();
      x2 = arrowsys_fast(d, c, b, alpha, y);
      t2.stop();
      t3.start();
      x3 = arrowsys_sparse<SparseLU<SparseMatrix<double>>>(d, c, b, alpha, y);
      t3.stop();
      t4.start();
      x4 = arrowsys_sparse<BiCGSTAB<SparseMatrix<double>>>(d, c, b, alpha, y);
      t4.stop();
    }
    times(i, 0) = n[i];
    times(i, 1) = t1.min();
    times(i, 2) = t2.min();
    times(i, 3) = t3.min();
    times(i, 4) = t4.min();
    times(i, 5) = (x4 - x3).norm();
  }
  return times;
}
/* SAM_LISTING_END_0 */


} // namespace arrowsys