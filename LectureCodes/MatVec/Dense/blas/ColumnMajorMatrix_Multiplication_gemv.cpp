/* ColumnMajorMatrix_Multiplication using the GEMV-routine from BLAS (+1 loop)
 */
#include "ColumnMajorMatrix.hpp"
ColumnMajorMatrix ColumnMajorMatrix::gemvMultiply(ColumnMajorMatrix &B) {
  assert(m == B.n);
  ColumnMajorMatrix C(n, B.m);  // important: must be zero: (done in
                                // constructor)
  double alpha(1.0), beta(1.0);
  for (int j = 0; j < m; ++j)
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, alpha, data, n, &(B(0, j)),
                1, beta, &C(0, j), 1);
  return C;
}
