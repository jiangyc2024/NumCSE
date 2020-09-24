/* ColumnMajorMatrix_Multiplication using the GEMM-routine from BLAS */
#include "ColumnMajorMatrix.h"
ColumnMajorMatrix ColumnMajorMatrix::gemmMultiply( ColumnMajorMatrix &B)
{
  assert(m==B.n);
  ColumnMajorMatrix C(n,B.m);//important: must be zero: (done in constructor)
  double alpha(1.0),beta(1.0);
  cblas_dgemm ( CblasColMajor, CblasNoTrans, CblasNoTrans, n, m, B.m, alpha, data, n, B.data, B.n, beta, C.data, C.n );
  return C;
}
