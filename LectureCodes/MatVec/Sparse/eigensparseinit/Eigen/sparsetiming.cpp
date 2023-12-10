///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s):
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "timetable.hpp"
#include <Eigen/Sparse>

using Eigen::SparseMatrix;
using Eigen::RowMajor;
using Eigen::Triplet;
using std::max;
using std::min;
using std::vector;

int main() {
  
  // Tabulate measured runtimes
  const vector<Eigen::Index> sizes {10, 100, 200, 800, 1000, 1500, 2000};

  SparseMatrix<double, RowMajor> spMat;
  const Eigen::Index bw{2};  // Half the bandwidth of band matrices
  Eigen::Index N = 0;

  /**
  * Using triplet list initialization.
  **/
  auto triplet_init = [&](){

    vector<Triplet<double>> triplets;
    triplets.reserve(N * 2 * bw);

    for (Eigen::Index i = 0; i < N; i++) {
      for (Eigen::Index j = max(0L, i - bw); j <= min(N - 1, i + bw); j++) {
        triplets.emplace_back(i, j, static_cast<double>(i + 1) * 0.5 + static_cast<double>(j + 1) * 0.7);
      }
    }

    spMat = SparseMatrix<double, RowMajor>(N, N);

    // The resulting matrix is already compressed
    spMat.setFromTriplets(triplets.begin(), triplets.end());
  };

  /**
   * Using mat.coeffRef(i, j) method NxN sparse matrix.
   * Reserve exactly as much as space as needed for each row.
   **/
  auto coefInsert_init1 = [&](){

    spMat = SparseMatrix<double, RowMajor>(N, N);
    spMat.reserve(vector<int>(N, 2 * bw + 1));

    for (Eigen::Index i = 0; i < N; i++) {
      for (Eigen::Index j = max(0L, i - bw); j <= min(N - 1, i + bw); j++) {
        spMat.coeffRef(i, j) += static_cast<double>(i + 1) * 0.5 + static_cast<double>(j + 1) * 0.7;
      }
    }

    spMat.makeCompressed();
  };

  /**
   * Using mat.coeffRef(i, j).
   * Reserve one element less for each row.
   * This is 100 times slower as compared to the previous version.
   * Takes ~3mins for 10^5 x 10^5 matrices.
   **/
  auto coefInsert_init2 = [&]() {

    spMat = SparseMatrix<double, RowMajor>(N, N);
    // Reserve one element less per row than actually needed.
    // In this case, it is really slow.
    spMat.reserve(std::vector<int>(N, 2 * bw));

    for (Eigen::Index i = 0; i < N; i++) {
      for (Eigen::Index j = max(0L, i - bw); j <= min(N - 1, i + bw); j++) {
        spMat.coeffRef(i, j) += static_cast<double>(i + 1) * 0.5 + static_cast<double>(j + 1) * 0.7;
      }
    }

    spMat.makeCompressed();
  };

  auto loopInit = [&](const Eigen::Index size) { 

    N = size; 
  };

  timeTable(sizes,
            {triplet_init, coefInsert_init1, coefInsert_init2}, loopInit,
            {"size", "triplets", "insert1", "insert2"});
}
