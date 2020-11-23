///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s):
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Sparse>
#include "timetable.hpp"

using namespace Eigen;
using namespace std;

SparseMatrix<double, RowMajor> spMat;
constexpr size_t bw{2};  // Half the bandwidth of band matrices
size_t N;

/**
 * Using mat.coeffRef(i, j) method NxN sparse matrix.
 * Reserve exactly as much as space as needed for each row.
 **/
void coefInsert_init1() {
  spMat = SparseMatrix<double, RowMajor>(N, N);
  spMat.reserve(vector<int>(N, 2 * bw + 1));

  for (size_t i = 0; i < N; i++)
    for (size_t j = max(0, (int)i - (int)bw); j <= min(N - 1, i + bw); j++)
      spMat.coeffRef(i, j) += (i + 1) * 0.5 + (j + 1) * 0.7;

  spMat.makeCompressed();
}

/**
 * Using mat.coeffRef(i, j).
 * Reserve one element less for each row.
 * This is 100 times slower as compared to the previous version.
 * Takes ~3mins for 10^5 x 10^5 matrices.
 **/
void coefInsert_init2() {
  spMat = SparseMatrix<double, RowMajor>(N, N);
  // Reserve one element less per row than actually needed.
  // In this case, it is really slow.
  spMat.reserve(std::vector<int>(N, 2 * bw));

  for (size_t i = 0; i < N; i++)
    for (size_t j = max(0, (int)i - (int)bw); j <= min(N - 1, i + bw); j++)
      spMat.coeffRef(i, j) += (i + 1) * 0.5 + (j + 1) * 0.7;

  spMat.makeCompressed();
}

/**
 * Using triplet list initialization.
 **/
void triplet_init() {
  vector<Triplet<double>> triplets;
  triplets.reserve(N * 2 * bw);

  for (size_t i = 0; i < N; i++)
    for (size_t j = max(0, (int)i - (int)bw); j <= min(N - 1, i + bw); j++)
      triplets.push_back({i, j, (i + 1) * 0.5 + (j + 1) * 0.7});

  spMat = SparseMatrix<double, RowMajor>(N, N);

  // The resulting matrix is already compressed
  spMat.setFromTriplets(triplets.begin(), triplets.end());
}

void loopInit(const size_t size) { N = size; }

int main() {
  // Tabulate measured runtimes
  vector<size_t> sizes {10, 100, 200, 800, 1000, 1500, 2000};
  timeTable( sizes,
            {triplet_init, coefInsert_init1, coefInsert_init2}, loopInit,
            {"size", "triplets", "insert1", "insert2"});
}
