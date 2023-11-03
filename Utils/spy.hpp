#pragma once

#include "matplotlibcpp.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>  //for nonZeros()

namespace plt = matplotlibcpp;

/**
 * @brief Equivalent to the old figure.h fig.spy(M) workflow
 * spy(M,title,fname) creates and saves a spyplot of M
 * @param title title of the plot
 * @param fname save as filename relative to working directory
 */
template <typename Matrix>
void spy(const Matrix &M, const std::string &title, const std::string &fname);

// Produce spy-plot of dense Eigen matrix
template <>
inline
void spy<Eigen::MatrixXd>(const Eigen::MatrixXd &M, const std::string &title,
                   const std::string &fname) {
  plt::figure();
  plt::spy(M, {{"marker", "o"}, {"markersize", "2"}, {"color", "b"}});
  plt::suptitle("nnz = " + std::to_string(
                               Eigen::SparseMatrix<double>(M.sparseView()).nonZeros()),
                {{"y", "0.09"}, {"fontsize", "small"}});
  plt::title(title);
  plt::savefig("./" + fname);
}

// Produce spy-plot of general Eigen matrix
template <typename Matrix>
void spy(const Matrix &M, const std::string &title, const std::string &fname) {
  spy(Eigen::MatrixXd(M), title, fname);
}