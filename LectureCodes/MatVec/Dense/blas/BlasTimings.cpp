#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <span>

#include "ColumnMajorMatrix.hpp"
#include "simpleTimer.hpp"

using std::cout;
using std::endl;

/* Main Routine for the timing of different
 * Matrix Matrix Multiplication implementations */
int main(int argc, char* const argv[]) {
  std::cout << "NumCSE timing code for BLAS routines" << std::endl;
  double T0(1e20);
  double T1(1e20);
  double T2(1e20);
  double T3(1e20);
  simpleTimer watch;
  int rep(5);
  int n(500);
  if (argc > 1) {
    auto args = std::span(argv, argc);
    n = static_cast<int>(strtol(args[1], nullptr, 10));
  }
  if (argc > 2) {
    auto args = std::span(argv, argc);
    rep = static_cast<int>(strtol(args[2], nullptr, 10));
  }
  // Declare Input Data
  ColumnMajorMatrix A(n, n);
  A.initRand();  // A.initGrow();
  ColumnMajorMatrix B(A);
  // The Results:
  ColumnMajorMatrix C(n, n);
  ColumnMajorMatrix D(n, n);
  ColumnMajorMatrix E(n, n);
  ColumnMajorMatrix F(n, n);
  // loop for repetitions (always take timing results over several
  // measurements!)
  for (int r = 0; r < rep; ++r) {
    watch.start();
    C = A.standardMultiply(B);
    T0 = std::min(T0, watch.getTime());
    watch.reset();
    watch.start();
    D = A.dotMultiply(B);
    T1 = std::min(T1, watch.getTime());
    watch.reset();
    watch.start();
    E = A.gemvMultiply(B);
    T2 = std::min(T2, watch.getTime());
    watch.reset();
    watch.start();
    F = A.gemmMultiply(B);
    T3 = std::min(T3, watch.getTime());
    watch.reset();
  }

  cout << "Timing Results: (min. of : " << rep << " Repetitions" << endl;
  cout << "N: " << n << " StraightForward: " << T0 << endl; 
  cout << "N: " << n << " dotMultiply: " << T1 << ", error: " << D.CalcErr(C) << endl; 
  cout << "N: " << n << " gemvMultiply: " << T2 << ", error: " << E.CalcErr(C) << endl; 
  cout << "N: " << n << " gemmMultiply: " << T3 << ", error: " << F.CalcErr(C) << endl; 
}
