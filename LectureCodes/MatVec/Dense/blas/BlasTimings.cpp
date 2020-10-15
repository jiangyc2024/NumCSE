#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "ColumnMajorMatrix.hpp"
#include "simpleTimer.hpp"
/* Main Routine for the timing of different
 * Matrix Matrix Multiplication implementations */
int main(int argc, char* const argv[]) {
  double T0(1e20), T1(1e20), T2(1e20), T3(1e20);
  simpleTimer watch;
  int rep(1), n(5);
  if (argc > 1) n = atoi(argv[1]);
  if (argc > 2) rep = atoi(argv[2]);
  // Declare Input Data
  ColumnMajorMatrix A(n, n);
  A.initRand();  // A.initGrow();
  ColumnMajorMatrix B(A);
  // The Results:
  ColumnMajorMatrix C(n, n), D(n, n), E(n, n), F(n, n);
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
  printf("Timing Results: (min. of : %i Repetitions) \n", rep);
  printf("N: %i StraightForward: %g \n", n, T0);
  printf("N: %i dotMultiply: %g ,error: %g \n", n, T1, D.CalcErr(C));
  printf("N: %i gemvMultiply: %g, error: %g \n", n, T2, E.CalcErr(C));
  printf("N: %i gemmMultiply: %g, error: %g \n", n, T3, F.CalcErr(C));
}
