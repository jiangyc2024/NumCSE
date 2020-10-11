#define daxpy_ cblas_daxpy
#include <iostream>

// Definition of the required BLAS function. This is usually done
// in a header file like \texttt{blas.h} that is included in the \eigen{}3
// distribution
extern "C" {
int daxpy_(const int* n, const double* da, const double* dx, const int* incx,
           double* dy, const int* incy);
}

using namespace std;

int main() {
  const int n = 5;     // length of vector
  const int incx = 1;  // stride
  const int incy = 1;  // stride
  double alpha = 2.5;  // scaling factor

  // Allocated raw arrays of doubles
  double* x = new double[n];
  double* y = new double[n];

  for (size_t i = 0; i < n; i++) {
    x[i] = 3.1415 * i;
    y[i] = 1.0 / (double)(i + 1);
  }

  cout << "x=[";
  for (size_t i = 0; i < n; i++) cout << x[i] << ' ';
  cout << "]" << endl;
  cout << "y=[";
  for (size_t i = 0; i < n; i++) cout << y[i] << ' ';
  cout << "]" << endl;

  // Call the BLAS library function passing pointers to all arguments
  // (Necessary when calling FORTRAN routines from C
  daxpy_(&n, &alpha, x, &incx, y, &incy);

  cout << "y = " << alpha << " * x + y = [";
  for (int i = 0; i < n; i++) cout << y[i] << ' ';
  cout << "]" << endl;
  return (0);
}
