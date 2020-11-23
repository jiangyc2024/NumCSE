///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "spdiags.hpp"
#include "spy.hpp"
#include "timetable.hpp"

using namespace std;
using namespace Eigen;

// compressed row storage
SparseMatrix<double, RowMajor> ACrs;

// compressed column storage (default in Eigen)
SparseMatrix<double, ColMajor> ACcs;

double tmp;  // prevent optimization
RowVectorXd v1;
VectorXd v2;
size_t n;

int main() {
  /* SAM_LISTING_BEGIN_0 */
  // Lambda, generates a sparse n x n matrix
  auto sparseM = [](const size_t n) {
    RowVectorXd diag_el(3);
    diag_el << -1, 2, 5;
    VectorXi diag_no(3);
    diag_no << -n / 2, 0, n / 2;
    MatrixXd B = diag_el.replicate(n, 1);
    // Place the columns of B along the diagonals specified by diag_no
    return spdiags(B, diag_no, n, n);
  };

  // Lambda, initializes environment for measurement with parameter i
  auto init = [sparseM](const size_t i) {
    n = std::pow(2, i);
    ACrs = ACcs = sparseM(n);
    v1 = RowVectorXd::Zero(n);
    v2 = VectorXd::Zero(n);
  };

  auto post = [] {
    tmp += v1 * v2;
  };  // use result of v1 and v2 after timing to prevent optimization

  auto rowAccessCcs = [] {
    for (int j = 1; j <= 5; ++j) v1 += (RowVectorXd)ACcs.row(n / 2);
  };

  auto colAccessCcs = [] {
    for (int j = 1; j <= 5; ++j) v2 += (VectorXd)ACcs.col(n / 2);
  };

  auto rowAccessCrs = [] {
    for (int j = 1; j <= 5; ++j) v1 += (RowVectorXd)ACrs.row(n / 2);
  };

  auto colAccessCrs = [] {
    for (int j = 1; j <= 5; ++j) v2 += (VectorXd)ACrs.col(n / 2);
  };

  spy(sparseM(16), "Pattern of sparse matrix for n=16",
      "spdiagsmatspy_cpp.eps");

  // tabulate access times for sparse matrix in CCS and CRS
  timeTable(VectorXi::LinSpaced(20, 1, 20),
            {rowAccessCcs, colAccessCcs, rowAccessCrs, colAccessCrs}, init,
            post, {"log(size)", "row CCS", "col CCS", "row CRS", "col CRS"});

  // plotting ...
  /* SAM_LISTING_END_0 */
  cout << tmp << endl;
}
