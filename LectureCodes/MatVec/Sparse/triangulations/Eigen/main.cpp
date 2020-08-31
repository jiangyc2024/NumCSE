///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <iostream>
#include "matplotlibcpp.h"
#include "pltextensions.hpp"
#include "spy.hpp"
#include "triangulation.hpp"

namespace plt = matplotlibcpp;

using namespace Eigen;

int main() {
  /* SAM_LISTING_BEGIN_0 */
  // Demonstration for visualizing a plane triangular mesh
  // Initialize node coordinates
  VectorXd x(10), y(10);
  // x and y coordinate of mesh
  x << 1.0, 0.60, 0.12, 0.81, 0.63, 0.09, 0.27, 0.54, 0.95, 0.96;
  y << 0.15, 0.97, 0.95, 0.48, 0.80, 0.14, 0.42, 0.91, 0.79, 0.95;
  // Then specify triangles through the indices of their vertices. These
  // indices refer to the ordering of the coordinates as given in the
  // vectors x and y.
  MatrixXi T(11, 3);
  T << 7, 1, 2, 5, 6, 2, 4, 1, 7, 6, 7, 2, 6, 4, 7, 6, 5, 0, 3, 6, 0, 8, 4, 3,
      3, 4, 6, 8, 1, 4, 9, 1, 8;
  plt::figure();
  plt::triplot(
      x, y, T,
      {{"color", "b"}, {"lw", ".5"}});  // drawing triangulation with numbers
  for (size_t i(0); i < x.size(); ++i) {
    plt::text(x[i], y[i], std::to_string(i), {{"ha", "center"}});
  }
  plt::plot(
      x, y,
      {{"marker", "*"}, {"color", "r"}, {"linestyle", " "}});  // mark vertices
  plt::savefig("./meshplot_cpp.eps");
  /* SAM_LISTING_END_0 */
  VectorXd x_start = x;
  VectorXd y_start = y;
  MatrixXi T_start = T;
  // Start new block
  {
    /* SAM_LISTING_BEGIN_1 */
    VectorXd x_ref, y_ref, xs, ys;
    MatrixXi T_ref;
    int refine_steps = 3;
    for (int i = 1; i <= refine_steps; ++i) {
      refinemesh(x, y, T, x_ref, y_ref, T_ref);
      plt::figure();
      plt::title("Refined mesh level " + std::to_string(i));
      plt::triplot(x_ref, y_ref, T_ref, {{"color", "b"}});
      plt::savefig("./rmesh" + std::to_string(i) + "_cpp.eps");
      smoothmesh(x_ref, y_ref, T_ref, xs, ys);
      x = xs;
      y = ys;
      T = T_ref;
      plt::figure();
      plt::title("Smoothed mesh level " + std::to_string(i));
      plt::triplot(x, y, T, {{"color", "b"}});
      plt::savefig("./smesh" + std::to_string(i) + "_cpp.eps");
    }
    /* SAM_LISTING_END_1 */
  }
  x = x_start;
  y = y_start;
  T = T_start;
  // Start new block
  {
    /* SAM_LISTING_BEGIN_2 */
    // Spyplot of graph laplacian
    VectorXd x_ref, y_ref, xs, ys;
    MatrixXi T_ref;
    int refine_steps = 3;
    for (int i = 1; i <= refine_steps; ++i) {
      refinemesh(x, y, T, x_ref, y_ref, T_ref);
      SparseMatrix<double> A_int;
      MatrixXd rhs;
      smoothmesh_analysis(x_ref, y_ref, T_ref, xs, ys, A_int, rhs);
      x = xs;
      y = ys;
      T = T_ref;

      MatrixXi E, Eb;
      // Extract the edge information of a mesh
      // E and Eb are matrices whose rows contain the numbers of the
      // endpoints of edges
      processmesh(T, E, Eb);
      int Nv = T.maxCoeff() + 1;
      int Ne = E.rows();
      int Nt = T.rows();
      std::string title = "level " + std::to_string(i) + ": " +
                          std::to_string(Nv) + " points, " +
                          std::to_string(Ne) + " edges, " + std::to_string(Nt) +
                          " cells";

      spy(A_int, title, "Lmat" + std::to_string(i) + "_cpp.eps");
    }
    /* SAM_LISTING_END_2 */
  }
  return 0;
}
