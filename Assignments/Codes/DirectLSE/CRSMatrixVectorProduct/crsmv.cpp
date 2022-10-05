/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: October 2021
 */

#include <Eigen/Dense>
#include <cassert>
#include <vector>
#include <iostream>

struct CRSMatrix {
  unsigned int m;                     // number of rows
  unsigned int n;                     // number of columns
  std::vector<double> val;            // value array
  std::vector<unsigned int> col_ind;  // same length as value array
  std::vector<unsigned int> row_ptr;  // length m+1, row_ptr[m] == val.size()
};

/* SAM_LISTING_BEGIN_1 */
template <typename VECTORTYPE_I, typename VECTORTYPE_II>
VECTORTYPE_I crsmv(const CRSMatrix &M, const VECTORTYPE_II &x) {
  assert((x.size() == M.n) && "Size mismatch between x and M");
  VECTORTYPE_I y(M.m);
  for (int k = 0; k < M.m; ++k) {
    y[k] = 0;
    for (int j = M.row_ptr[k]; j < M.row_ptr[k + 1]; ++j) {
      y[k] += M.val[j] * x[M.col_ind[j]];
    }
  }
  return y;
}
/* SAM_LISTING_END_1 */

int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "C++ program for course Numerical Methods for CSE" << std::endl;
  std::cout << "Multiplication of CRS matrix with vector" << std::endl;

  CRSMatrix M;
  M.m = 6;
  M.n = 6;
  M.val = {10.0, -2.0, 3.0, 9.0, 3.0, 7.0,  8.0, 7.0, 3.0, 8.0,
           7.0,  5.0,  8.0, 9.0, 9.0, 13.0, 4.0, 2.0, -1.0};
  M.col_ind = {0, 4, 0, 1, 5, 1, 2, 3, 0, 2, 3, 4, 1, 3, 4, 5, 1, 4, 5};
  M.row_ptr = {0, 2, 5, 8, 12, 16, 19};

  Eigen::MatrixXd A(6, 6);
  // clang format off 
  A << 10, 0, 0, 0, -2, 0,
       3, 9, 0, 0, 0, 3,
       0, 7, 8, 7, 0, 0,
       3, 0, 8, 7, 5, 0,
       0, 8, 0, 9, 9, 13,
       0, 4, 0, 0, 2, -1;
  // clang format on 
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(6,1.0,6.0);
  std::cout << crsmv<Eigen::VectorXd>(M,x).transpose() << std::endl;
  std::cout << (A*x).transpose() << std::endl;
  return 0;
}
