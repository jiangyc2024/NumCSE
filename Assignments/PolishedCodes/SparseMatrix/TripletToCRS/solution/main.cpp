#include <iostream>

#include "timer.h"
#include "triplettoCRS.hpp"

int main() {
  // Initialization
  constexpr std::size_t nrows = 7, ncols = 5, ntriplets = 9;
  TripletMatrix<double> T;
  CRSMatrix<double> C, D;

  // TODO: (2-13.f) Construct T here.
  // START
  T.rows = nrows;
  T.cols = ncols;
  T.triplets.reserve(ntriplets);  // Always reserve space if you can
  // END

  for (std::size_t i = 0; i < ntriplets; ++i) {
    // TODO: (2-13.f) Use this loop to push back triplets in your matrix. Insert
    // triplet with arguments: (rand() % nrows, rand() % ncols, rand() % 1000))
    // START
    // Test unordered triplets, random and possibly repeated triplets
    T.triplets.push_back(
        Triplet<double>(rand() % nrows, rand() % ncols, rand() % 1000));
    // END
  }

  std::cout << "***Test conversion with random matrices***" << std::endl;
  tripletToCRS_insertsort(T, C);
  tripletToCRS_sortafter(T, D);
  // TODO: (2-13.f) If you implemented densify(), compute the Frobenius norm of
  // $T - C$.
  // START
  std::cout << "--> Frobenius norm of T - C: "
            << (T.densify() - C.densify()).norm() << std::endl;
  std::cout << "--> Frobenius norm of T - D: "
            << (T.densify() - D.densify()).norm() << std::endl;
  // END

  std::cout << "Enter 0 for no runtime measurements." << std::endl;
  int noruntime = 0;
  std::cin >> noruntime;
  if (!noruntime) {
    return 0;
  }

  std::cout << "***Runtime test***" << std::endl;

  // Play around with parameters (also introducing lambda functions)
  auto frows = [](std::size_t M) { return 2 * M; };
  auto fcols = [](std::size_t M) { return M; };
  auto ftriplets = [](std::size_t M) { return M * 5; };

  // Runtime test
  Timer insertsort_timer, sortafter_timer;
  for (std::size_t n = 2; n < 1024; n *= 2) {
    std::cout << "Runtime for " << n << "x" << n / 2
              << " matrix (with nnz(A) <= " << n << "):" << std::endl;
    TripletMatrix<double> A;
    CRSMatrix<double> E, F;

    // TODO: (2-13.f) Construct and initialize A for runtime measurements.
    // START
    A.rows = frows(n);  // nrows
    A.cols = fcols(n);  // ncols
    A.triplets.reserve(ftriplets(n));
    for (std::size_t i = 0; i < ftriplets(n); ++i) {
      A.triplets.push_back(
          Triplet<double>(rand() % A.rows, rand() % A.cols, rand() % 1000));
    }
    // END

    insertsort_timer.start();
    tripletToCRS_insertsort(A, E);
    insertsort_timer.stop();

    sortafter_timer.start();
    tripletToCRS_sortafter(A, F);
    sortafter_timer.stop();
    std::cout << "InsertSort took: " << insertsort_timer.duration() << " s."
              << std::endl;
    std::cout << "SortAfter took:  " << sortafter_timer.duration() << " s."
              << std::endl;
  }
}