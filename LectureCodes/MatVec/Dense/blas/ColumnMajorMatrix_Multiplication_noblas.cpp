/* ColumnMajorMatrix_Multiplication: straightforward standard implementation for
 * a Matrix Multiplication (3 loops) */
#include "ColumnMajorMatrix.hpp"
ColumnMajorMatrix ColumnMajorMatrix::standardMultiply(ColumnMajorMatrix &B) {
  assert(m == B.n);             // only support square matrices
  ColumnMajorMatrix C(n, B.m);  // important: must be zero: (done in
                                // constructor)
  for (int j = 0; j < B.m; ++j)
    for (int i = 0; i < n; ++i)
      for (int k = 0; k < m; ++k) C(i, j) += operator()(i, k) * B(k, j);
  return C;
}
