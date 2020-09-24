/* ColumnMajorMatrix_Multiplication.cpp using the DOT-routine from BLAS (and 2 loops) */
#include "ColumnMajorMatrix.h"
ColumnMajorMatrix ColumnMajorMatrix::dotMultiply( ColumnMajorMatrix &B)
{
  assert(m==B.n); 
  ColumnMajorMatrix C(n,B.m);
  for (int j=0;j<B.m;++j)
    for(int i=0;i<n;++i)
      C(i,j)=cblas_ddot(this->m, &(operator()(i,0)) , n, &(B(0,j)), 1);
  return C;
}
