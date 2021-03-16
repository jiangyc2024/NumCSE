/*  Author: Manfred Quack
 *  ColumnMajorMatrix.h
 *  This Class Implements a ColumnMajor Matrix Structure in C++
 *  - it provides an access operator () using ColumnMajor Access
 *  - it also provides 4 different implementations of a Matrix-Matrix
 * Multiplication
 */
#include <stdio.h>

#include <cassert>
#include <cstdlib>
#include <iostream>
#ifndef _USE_MKL
#ifdef _MAC_OS
#include <Accelerate/Accelerate.h>
#endif
#ifdef _LINUX
extern "C" {
#include <cblas.h>
}
#endif
#else
#include <mkl.h>
#endif
#include <cmath>
typedef double Real;
using namespace std;

class ColumnMajorMatrix {
 private:  // Data Members:
  Real *data;
  int n, m;

 public:
  //------Class Managment---------:
  // Constructors
  ColumnMajorMatrix(int _n, int _m);
  ColumnMajorMatrix(const ColumnMajorMatrix &B);
  // Destructor
  ~ColumnMajorMatrix();
  // Assignment operator
  ColumnMajorMatrix &operator=(const ColumnMajorMatrix &B);
  // Access for ColumnMajorArrangement
  inline Real &operator()(int i, int j) {
    assert(i < n && j < m);
    return data[n * j + i];
  }
  //----Different Implementations for the Multiplication----/
  // All of these implementations have only been checked for  Square Matrices
  // straigthforward implementation for a Matrix Multiplication
  ColumnMajorMatrix standardMultiply(ColumnMajorMatrix &B);
  // Implementations using DOT-, GEMV- and GEMM from BLAS:
  ColumnMajorMatrix dotMultiply(ColumnMajorMatrix &B);
  ColumnMajorMatrix gemvMultiply(ColumnMajorMatrix &B);
  ColumnMajorMatrix gemmMultiply(ColumnMajorMatrix &B);
  //------Other Functions---------:
  // Function to initialize matrix with 1,2,3...
  void initGrow();
  void initRand();
  void print();
  Real CalcErr(const ColumnMajorMatrix &B);
};
