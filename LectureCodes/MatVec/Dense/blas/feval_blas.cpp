#include "feval_blas.h"
// same as in feval_noblas but now with BLAS-operations
feval_blas::feval_blas(int _n, const double* _A, const double* _B,
                       const double* _u, const double* _v)
    : n(_n), M(_n, _n) {
  const double alpha = 1.0;
  const double beta = 1.0;
  double* temp_v1 = new double[n];
  double* temp_v2 = new double[n];
  for (int i = 0; i < n; ++i) {
    temp_v1[i] = 0;
    temp_v2[i] = 0;
  }
  // calculate tempv1=B*u
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, alpha, _B, n, _u, 1, beta,
              temp_v1, 1);
  // calculate tempv2=A^T*v (which is the same as uv^TA)
  cblas_dgemv(CblasColMajor, CblasTrans, n, n, alpha, _A, n, _v, 1, beta,
              temp_v2, 1);
  // calculate M+=u*temp_v2 (which is: M+=u*v^TA)
  // BLAS: DGER - perform the rank 1 operation   A := alpha*x*y' + A
  cblas_dger(CblasColMajor, n, n, alpha, _u, 1, temp_v2, 1, &(M(0, 0)), n);
  // calculate M+=v*temp_v1 (which is: M+=B*u*v^T)
  cblas_dger(CblasColMajor, n, n, alpha, temp_v1, 1, _v, 1, &(M(0, 0)), n);

  delete[] temp_v1;
  delete[] temp_v2;
}

void feval_blas::eval(const double* x, double* y) const {
  assert(x != NULL);
  assert(y != NULL);
  M.blasVectorMultiply(x, y);
}
