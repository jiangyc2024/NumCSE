///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s):
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Sparse>
#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std::chrono;

/**
 * Using mat.coeffRef(i, j) method NxN sparse matrix.
 * Reserve exactly as much as space as needed for each row.
 **/
std::size_t coefInsert_init1(const std::size_t N, const std::size_t bw) {
  auto tic = high_resolution_clock::now();

  SparseMatrix<double, RowMajor> spMat(N, N);
  spMat.reserve(std::vector<int>(N, 2 * bw + 1));

  for (std::size_t i = 0; i < N; i++)
    for (std::size_t j = std::max(0, (int)i - (int)bw);
         j <= std::min(N - 1, i + bw); j++)
      spMat.coeffRef(i, j) += (i + 1) * 0.5 + (j + 1) * 0.7;

  spMat.makeCompressed();

  auto toc = high_resolution_clock::now();
  return duration_cast<milliseconds>(toc - tic).count();
}

/**
 * Using mat.coeffRef(i, j).
 * Reserve one element less for each row.
 * This is 100 times slower as compared to the previous version.
 * Takes ~3mins for 10^5 x 10^5 matrices.
 **/
std::size_t coefInsert_init2(const std::size_t N, const std::size_t bw) {
  auto tic = high_resolution_clock::now();
  SparseMatrix<double, RowMajor> spMat(N, N);

  // Reserve one element less per row than actually needed.
  // In this case, it is really slowly.
  spMat.reserve(std::vector<int>(N, 2 * bw));

  for (std::size_t i = 0; i < N; i++)
    for (std::size_t j = std::max(0, (int)i - (int)bw);
         j <= std::min(N - 1, i + bw); j++)
      spMat.coeffRef(i, j) += (i + 1) * 0.5 + (j + 1) * 0.7;

  spMat.makeCompressed();

  auto toc = high_resolution_clock::now();
  return duration_cast<milliseconds>(toc - tic).count();
}

/**
 * Using triplet list initialization.
 **/
std::size_t triplet_init(const std::size_t N, const std::size_t bw) {
  auto tic = high_resolution_clock::now();

  std::vector<Triplet<double>> triplets;
  triplets.reserve(N * 2 * bw);

  for (std::size_t i = 0; i < N; i++)
    for (std::size_t j = std::max(0, (int)i - (int)bw);
         j <= std::min(N - 1, i + bw); j++)
      // triplets.push_back(Triplet<double>(i, j, (i+1)*0.5+(j+1)*0.7));
      triplets.push_back({i, j, (i + 1) * 0.5 + (j + 1) * 0.7});

  SparseMatrix<double, RowMajor> spMat(N, N);
  spMat.setFromTriplets(triplets.begin(), triplets.end());

  auto toc = high_resolution_clock::now();
  return duration_cast<milliseconds>(toc - tic).count();
}

int main() {
  // Half the bandwidth of band matrices
  std::size_t bw = 2;
  std::cout << std::scientific << std::setprecision(3);
  for (std::size_t n = 1; n <= 5; n++) {
    std::size_t N = std::pow(10, n);
    std::size_t time_triplet = triplet_init(N, bw);
    std::size_t time_insert1 = coefInsert_init1(N, bw);
    // std::size_t time_insert2 = coefInsert_init2(N, bw);

    std::cout << n << "," << time_triplet << "," << time_insert1 << ","
              << std::endl;
    // << time_insert2 << std::endl;
  }
  return 0;
}
