///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <algorithm>
#include <iostream>
#include <list>
#include <numeric>

#include <Eigen/Dense>
#include <Eigen/Sparse>

inline
/* SAM_LISTING_BEGIN_0 */
//! @brief \eigen function extracting the edge information of a mesh
//! @param[in] T Matrix $M \times 3$ containing the vertex numbers of each of
//! the $M$ triangles of a triangular mesh
//! @param[out] E Matrix $\sharp \mathcal{E} \times 2$ whose rows correspond to
//! the edges of the mesh and give the indices (corresponding to T) of the
//! endpoints of the edges
//! @param[out] Eb Matrix $\sharp \mathcal{\Gamma} \times 2$ whose rows
//! correspond to the boundary edges of the mesh and give the indices
//! (corresponding to T) of the endpoints of the edges
void processmesh(const Eigen::MatrixXi& T, Eigen::MatrixXi& E,
                 Eigen::MatrixXi& Eb) {
  // Number of nodes of the triangular mesh
  const int N = T.maxCoeff() + 1;
  // Number of triangles of the mesh
  const Eigen::Index M = T.rows();
  // Triplet vector for initializing a sparse matrix
  std::vector<Eigen::Triplet<int> > triplets;
  // Reserve enough space for the vector to prevent reallocation
  triplets.reserve(3 * M);
  // Loop over all triangles to get all edges at least once
  for (int k = 0; k < M; ++k) {
    /* Loop over all possible  edge combinations of a triangle and insert them
     *     2
     *    / \		--> Edges: [0,1], [0,2] [1,2]
     *   /   \
     *  0-----1
     */
    for (int i = 0; i < 2; ++i) {
      for (int j = i + 1; j < 3; ++j) {
        if (T(k, i) < T(k, j)) { 
          // insert combinations sorted wrt. to node numbering (ascending)
          triplets.emplace_back(T(k, i), T(k, j), 1);
        } 
        else {
          triplets.emplace_back(T(k, j), T(k, i), 1);
        }
      }
    }
  }
  // Initialize sparse matrix
  Eigen::SparseMatrix<int> A(N, N);
  // Build sparse Matrix A from Triplet vector
  // The rows and the columns of the matrix correspond to nodes of
  // the mesh. If A(i,j) is different from zero the two nodes with
  // numbers i and j are connected by an edge. Note that A(i,j) is
  // equal to 2, if [i,j] is an interior edge. For a boundary edge
  // [i,j] A(i,j) == 1, because an interior edge belongs only to
  // one triangle.
  A.setFromTriplets(triplets.begin(), triplets.end());
  // Get the number of nnz entries which is corresponding to the
  // number of edges
  const Eigen::Index E_size = A.nonZeros();
  // resize E accordingly
  E.resize(E_size, 2);
  // Used to dynamically store the boundary edges {startpoint1, endpoint1, ...}
  std::vector<int> eb;
  // Iteration over non-zero elements of the sparse matrix
  // works only for col major
  int E_counter = 0;
  for (int i = 0; i < A.outerSize(); ++i) {
    for (Eigen::SparseMatrix<int>::InnerIterator it(A, i); it; ++it) {
      if (it.value() == 1) {  // boundary edge == 1
        eb.push_back(static_cast<int>(it.row()));
        eb.push_back(static_cast<int>(it.col()));
      }
      // boundary edge == 1 or interior edge == 2
      E(E_counter, 0) = static_cast<int>(it.row());
      E(E_counter, 1) = static_cast<int>(it.col());
      ++E_counter;
    }
  }
  // Use the \eigen Map function to assign the created vector to the
  // matrices Eb. Note that the RowMajor keyword. Copying is
  // unfortunately not prevented
  Eb = Eigen::Matrix<int, -1, -1, Eigen::RowMajor>::Map(eb.data(),
         static_cast<int>(eb.size()) / 2, 2);
}
/* SAM_LISTING_END_0 */

inline
/* SAM_LISTING_BEGIN_1 */
//! @brief \eigen function extracting the triangle-edge mapping
//! @param[in] T Matrix $M \times 3$ containing the vertex numbers of each of
//! the M triangles of a triangular mesh
//! @param[in] E Matrix $\sharp \mathcal{E} \times 2$ whose rows correspond to
//! the edges of the mesh and give the indices (corresponding to T) of the
//! endpoints of the edges
//! @param[out] ET Matrix $N \times 3$ whose rows contain the index numbers of
//! the edges of each triangle
void getinfo(const Eigen::MatrixXi& T, const Eigen::MatrixXi& E,
             Eigen::MatrixXi& ET) {
  // Number of edges $ L = \sharp \mathcal{E} $
  const Eigen::Index L = E.rows();
  // Number of nodes of the triangular mesh
  const int N = T.maxCoeff() + 1;
  // Triplet vector for initializing a sparse matrix
  std::vector<Eigen::Triplet<int> > triplets;
  // Reserve enough space for the vector to prevent reallocation
  triplets.reserve(2 * L);
  // loop over all edges
  for (int i = 0; i < L; ++i) {
    triplets.emplace_back(E(i, 0), E(i, 1), i);
    triplets.emplace_back(E(i, 1), E(i, 0), i);  // symmetrical
  }
  // Initialize sparse matrix
  Eigen::SparseMatrix<int> A(N, N);
  // Build sparse Matrix A from Triplet vector
  // Create a sparse $ N \times N $ A matrix, for which A(i,j) is zero
  // if the nodes i and j are not connected by an edge, and for which
  // A(i,j) gives the index of the edge from i to j, if it exists.
  // The index of an edge corresponds to its row number in the
  // E matrix.
  A.setFromTriplets(triplets.begin(), triplets.end());
  // resize ET accordingly
  ET.resize(T.rows(), 3);
  // Get the 3 edges corresponding to a triangle
  for (int i = 0; i < T.rows(); ++i) {
    ET(i, 0) = A.coeff(T(i, 1), T(i, 2));
    ET(i, 1) = A.coeff(T(i, 0), T(i, 2));
    ET(i, 2) = A.coeff(T(i, 0), T(i, 1));
  }
}
/* SAM_LISTING_END_1 */

inline
/* SAM_LISTING_BEGIN_2 */
//! @brief \eigen function extracting the triangle edge mapping
//! @param[in] x Vector of dim $N$ containing the x-coordinates of nodes of the
//! mesh
//! @param[in] y Vector of dim $N$ containing the y-coordinates of nodes of the
//! mesh
//! @param[in] T Matrix $M \times 3$ containing the vertex numbers of each of
//! the M triangles of a triangular mesh
//! @param[out] x_ref Vector of dim $N$ containing the x-coordinates of nodes of
//! the refined mesh
//! @param[out] y_ref Vector of dim $N$ containing the y-coordinates of nodes of
//! the refined mesh
//! @param[out] T_ref Matrix $M \times 3$ containing the vertex numbers of each
//! of the M triangles of a triangular refined mesh
void refinemesh(const Eigen::VectorXd& x, const Eigen::VectorXd& y,
                const Eigen::MatrixXi& T, Eigen::VectorXd& x_ref,
                Eigen::VectorXd& y_ref, Eigen::MatrixXi& T_ref) {
  // Instantiate needed matrices for passing by reference
  Eigen::MatrixXi E; 
  Eigen::MatrixXi Eb;
  Eigen::MatrixXi ET;
  // Extract the edge information of a mesh
  // E and Eb are matrices whose rows contain the numbers of the
  // endpoints of edges
  processmesh(T, E, Eb);
  // resize the refined coordinates and fill them with the old ones
  x_ref.resize(x.size() + E.rows());
  y_ref.resize(y.size() + E.rows());
  x_ref.head(x.size()) = x;
  y_ref.head(y.size()) = y;
  // Add midpoints of edges to the coordinate vectors (refined coord)
  for (int i = 0; i < E.rows(); ++i) {
    x_ref(i + x.size()) = (x(E(i, 0)) + x(E(i, 1))) / 2.0;
    y_ref(i + y.size()) = (y(E(i, 0)) + y(E(i, 1))) / 2.0;
  }
  // Extract the triangle-edge mapping
  getinfo(T, E, ET);
  // Build a new list of triangles
  const Eigen::Index Nt = T.rows();  // Number of triangles
  const int Nv = static_cast<int>(x.size());  // Number of vertices
  // Resize the refined triangle (refinement creates 4 triangles per
  // 1 triangle
  T_ref.resize(4 * Nt, 3);
  // Fill T_ref with the new triangles
  for (Eigen::Index i = 0; i < Nt; ++i) {
    // 1st son triangle
    T_ref(4 * i, 0) = T(i, 0);        // old vertex
    T_ref(4 * i, 1) = ET(i, 1) + Nv;  // new vertex
    T_ref(4 * i, 2) = ET(i, 2) + Nv;  // new vertex
    // 2nd son triangle
    T_ref(4 * i + 1, 0) = T(i, 1);
    T_ref(4 * i + 1, 1) = ET(i, 0) + Nv;
    T_ref(4 * i + 1, 2) = ET(i, 2) + Nv;
    // 3rd son triangle
    T_ref(4 * i + 2, 0) = T(i, 2);
    T_ref(4 * i + 2, 1) = ET(i, 0) + Nv;
    T_ref(4 * i + 2, 2) = ET(i, 1) + Nv;
    // 3th (interior) son triangle
    T_ref(4 * i + 3, 0) = ET(i, 0) + Nv;
    T_ref(4 * i + 3, 1) = ET(i, 1) + Nv;
    T_ref(4 * i + 3, 2) = ET(i, 2) + Nv;
  }
}
/* SAM_LISTING_END_2 */

inline
/* SAM_LISTING_BEGIN_3 */
//! @brief \eigen function for smoothing a planar triangular mesh by moving all
//! vertices to the barycenter of their neighboring nodes.
//! @param[in] x Vector of dim $N$ containing the x-coordinates of nodes of the
//! mesh
//! @param[in] y Vector of dim $N$ containing the y-coordinates of nodes of the
//! mesh
//! @param[in] T Matrix $M \times 3$ containing the vertex numbers of each of
//! the M triangles of a triangular mesh
//! @param[out] xs Vector of dim $N$ containing the x-coordinates at the
//! barycenter of their neighboring nodes
//! @param[out] ys Vector of dim $N$ containing the y-coordinates at the
//! barycenter of their neighboring nodes
void smoothmesh(const Eigen::VectorXd& x, const Eigen::VectorXd& y,
                const Eigen::MatrixXi& T, Eigen::VectorXd& xs,
                Eigen::VectorXd& ys) {
  // Number of nodes of the mesh
  const Eigen::Index Nv = x.size();
  // Instantiate needed matrices for passing by reference
  Eigen::MatrixXi E;
  Eigen::MatrixXi Eb;
  // Extract the edge information of a mesh
  // E and Eb are matrices whose rows contain the numbers of the
  // endpoints of edges
  processmesh(T, E, Eb);
  // Initialize a vector containing all boundary nodes indices
  std::vector<int> bd_nodes(Eb.reshaped().begin(), Eb.reshaped().end());
  // erase the duplicates
  std::sort(bd_nodes.begin(), bd_nodes.end());
  auto last = std::unique(bd_nodes.begin(), bd_nodes.end());
  bd_nodes.erase(last, bd_nodes.end());
  // Number of boundary nodes
  const Eigen::Index Nb = bd_nodes.size();
  // get the coordinates of the boundary nodes
  Eigen::VectorXd x_bd(Nb);
  Eigen::VectorXd y_bd(Nb);
  int counter = 0;
  for (auto idx : bd_nodes) {
    x_bd(counter) = x(idx);
    y_bd(counter) = y(idx);
    ++counter;
  }
  // Get interior nodes
  std::vector<int> int_nodes;
  // Number of interior nodes
  const Eigen::Index Ni = Nv - Nb;
  int_nodes.reserve(Ni);
  // Create a vector with {$0,1,2, \ldots, \texttt{Nv}-1$}
  std::vector<int> nodes_map(Nv);
  std::iota(nodes_map.begin(), nodes_map.end(), 0);
  // Use the function set_difference to obtain the interior nodes
  std::set_difference(nodes_map.begin(), nodes_map.end(), bd_nodes.begin(),
                      bd_nodes.end(),
                      std::inserter(int_nodes, int_nodes.begin()));
  // Inverse map from new positions to old position/index in the matrix
  // nodes_map_inv[i] returns old position/index of i in the matrix, wheras
  // i is the new position/index
  // nodes_map is partioned in [interior indices, extorior indices]
  std::vector<int> nodes_map_inv(int_nodes.begin(), int_nodes.end());
  nodes_map_inv.insert(nodes_map_inv.end(), bd_nodes.begin(), bd_nodes.end());
  // Map to original/old node position
  // reuse the above vector nodes_map
  // a lambda function is used in the sort function to get the map
  std::sort(nodes_map.begin(), nodes_map.end(),
            [&nodes_map_inv](size_t a, size_t b) {
              return nodes_map_inv[a] < nodes_map_inv[b];
            });
  // Triplet vectors for initializing a sparse matrix
  std::vector<Eigen::Triplet<int> > triplets_int;  // interior
  std::vector<Eigen::Triplet<int> > triplets_bd;   // exterior
  // Reserve space for the interior nodes following from the
  // euler characteristics, overestimation since boundary of
  // interior triangulation is not known
  triplets_int.reserve(2 * (3 * Ni - 3));
  // Reserve some space for boundary part
  triplets_bd.reserve(Nb);
  // Keeps track of $\sharp S(i)$
  Eigen::VectorXi neighbours(Nv);
  neighbours.setZero();
  // Assemble the (symmetric) interior graph laplacian (triplets_int)
  // and the exterior graph laplacian (not symmetric)
  for (int i = 0; i < E.rows(); ++i) {
    const int idx1 = nodes_map[E(i, 0)];
    const int idx2 = nodes_map[E(i, 1)];
    if (idx1 < Ni && idx2 < Ni) {  // edge belongs to interior area
      triplets_int.emplace_back(idx1, idx2, -1);
      triplets_int.emplace_back(idx2, idx1, -1);
    } else if (idx1 < Ni) {  // only 1st edge node belongs to interior area
      triplets_bd.emplace_back(idx1, idx2 - Ni, -1);
    } else if (idx2 < Ni) {  // only 2nd edge node belongs to interior area
      triplets_bd.emplace_back(idx2, idx1 - Ni, -1);
    }
    ++neighbours(idx1);
    ++neighbours(idx2);
  }
  // Insert $\sharp S(i)$ on the diagonal of the interior matrix
  for (int i = 0; i < Ni; ++i) { // interior
    triplets_int.emplace_back(i, i, neighbours(i));
  }
  // Build matrices from Triplets
  Eigen::SparseMatrix<double> A_int(Ni, Ni);
  A_int.setFromTriplets(triplets_int.begin(), triplets_int.end());
  Eigen::SparseMatrix<double> A_bd(Ni, Nb);
  A_bd.setFromTriplets(triplets_bd.begin(), triplets_bd.end());
  // Compute rhs -A_bd*x and -A_bd*y and concatenate them to allow
  // solving only once
  Eigen::MatrixXd rhs(x_bd.size(), 2);
  rhs << x_bd, y_bd;
  rhs.applyOnTheLeft(-A_bd);
  // Instatiate solver
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
  // Indicate that the pattern of the input matrix is symmetric
  solver.isSymmetric(true);
  // Compute decomposition
  solver.compute(A_int);
  // Solve the linear system of equations with the corresponding rhs
  Eigen::MatrixXd xy_int = solver.solve(rhs);
  // Back transformation of boundary nodes (still in the same place)
  xs = x;
  ys = y;
  // Back transformation with the invers mapping
  for (int i = 0; i < Ni; ++i) {
    xs(nodes_map_inv[i]) = xy_int(i, 0);
    ys(nodes_map_inv[i]) = xy_int(i, 1);
  }
}
/* SAM_LISTING_END_3 */

inline
/* SAM_LISTING_BEGIN_4 */
//! @brief Copy of smoothmesh for the analysis of the runtimg, see ADDED FOR
//! ANALYIS
//! @param[out] additionally returns A_int and rhs for timings or spy-plots
void smoothmesh_analysis(const Eigen::VectorXd& x, const Eigen::VectorXd& y,
                         const Eigen::MatrixXi& T, Eigen::VectorXd& xs,
                         Eigen::VectorXd& ys,
                         Eigen::SparseMatrix<double>& A_int_out,
                         Eigen::MatrixXd& rhs_out) {
  // Number of nodes of the mesh
  const Eigen::Index Nv = x.size();
  // Instantiate needed matrices for passing by reference
  Eigen::MatrixXi E;
  Eigen::MatrixXi Eb;
  // Extract the edge information of a mesh
  // E and Eb are matrices whose rows contain the numbers of the
  // endpoints of edges
  processmesh(T, E, Eb);
  // Initialize a vector containing all boundary nodes indices
  std::vector<int> bd_nodes(Eb.reshaped().begin(), Eb.reshaped().end());
  // erase the duplicates
  std::sort(bd_nodes.begin(), bd_nodes.end());
  auto last = std::unique(bd_nodes.begin(), bd_nodes.end());
  bd_nodes.erase(last, bd_nodes.end());
  // Number of boundary nodes
  const Eigen::Index Nb = bd_nodes.size();
  // get the coordinates of the boundary nodes
  Eigen::VectorXd x_bd(Nb);
  Eigen::VectorXd y_bd(Nb);
  int counter = 0;
  for (auto idx : bd_nodes) {
    x_bd(counter) = x(idx);
    y_bd(counter) = y(idx);
    ++counter;
  }
  // Get interior nodes
  std::vector<int> int_nodes;
  // Number of interior nodes
  const Eigen::Index Ni = Nv - Nb;
  int_nodes.reserve(Ni);
  // Create a vector with {$0,1,2, \ldots, \texttt{Nv}-1$}
  std::vector<int> nodes_map(Nv);
  std::iota(nodes_map.begin(), nodes_map.end(), 0);
  // Use the function set_difference to obtain the interior nodes
  std::set_difference(nodes_map.begin(), nodes_map.end(), bd_nodes.begin(),
                      bd_nodes.end(),
                      std::inserter(int_nodes, int_nodes.begin()));
  // Inverse map from new positions to old position/index in the matrix
  // nodes_map_inv[i] returns old position/index of i in the matrix, wheras
  // i is the new position/index
  // nodes_map is partioned in [interior indices, extorior indices]
  std::vector<int> nodes_map_inv(int_nodes.begin(), int_nodes.end());
  nodes_map_inv.insert(nodes_map_inv.end(), bd_nodes.begin(), bd_nodes.end());
  // Map to original/old node position
  // reuse the above vector nodes_map
  // a lambda function is used in the sort function to get the map
  std::sort(nodes_map.begin(), nodes_map.end(),
            [&nodes_map_inv](size_t a, size_t b) {
              return nodes_map_inv[a] < nodes_map_inv[b];
            });
  // Triplet vectors for initializing a sparse matrix
  std::vector<Eigen::Triplet<int> > triplets_int;  // interior
  std::vector<Eigen::Triplet<int> > triplets_bd;   // exterior
  // Reserve space for the interior nodes following from the
  // euler characteristics, overestimation since boundary of
  // interior triangulation is not known
  triplets_int.reserve(2 * (3 * Ni - 3));
  // Reserve some space for boundary part
  triplets_bd.reserve(Nb);
  // Keeps track of $\sharp S(i)$
  Eigen::VectorXi neighbours(Nv);
  neighbours.setZero();
  // Assemble the (symmetric) interior graph laplacian (triplets_int)
  // and the exterior graph laplacian (not symmetric)
  for (int i = 0; i < E.rows(); ++i) {
    const int idx1 = nodes_map[E(i, 0)];
    const int idx2 = nodes_map[E(i, 1)];
    if (idx1 < Ni && idx2 < Ni) {  // edge belongs to interior area
      triplets_int.emplace_back(idx1, idx2, -1);
      triplets_int.emplace_back(idx2, idx1, -1);
    } else if (idx1 < Ni) {  // only 1st edge node belongs to interior area
      triplets_bd.emplace_back(idx1, idx2 - Ni, -1);
    } else if (idx2 < Ni) {  // only 2nd edge node belongs to interior area
      triplets_bd.emplace_back(idx2, idx1 - Ni, -1);
    }
    ++neighbours(idx1);
    ++neighbours(idx2);
  }
  // Insert $\sharp S(i)$ on the diagonal of the interior matrix
  for (int i = 0; i < Ni; ++i)  // interior
    triplets_int.emplace_back(i, i, neighbours(i));
  // Build matrices from Triplets
  Eigen::SparseMatrix<double> A_int(Ni, Ni);
  A_int.setFromTriplets(triplets_int.begin(), triplets_int.end());
  Eigen::SparseMatrix<double> A_bd(Ni, Nb);
  A_bd.setFromTriplets(triplets_bd.begin(), triplets_bd.end());
  // Compute rhs -A_bd*x and -A_bd*y and concatenate them to allow
  // solving only once
  Eigen::MatrixXd rhs(x_bd.size(), 2);
  rhs << x_bd, y_bd;
  rhs.applyOnTheLeft(-A_bd);
  // Instatiate solver
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
  // Indicate that the pattern of the input matrix is symmetric
  solver.isSymmetric(true);
  // Compute decomposition
  solver.compute(A_int);
  // Solve the linear system of equations with the corresponding rhs
  Eigen::MatrixXd xy_int = solver.solve(rhs);
  // Back transformation of boundary nodes (still in the same place)
  xs = x;
  ys = y;
  // Back transformation with the invers mapping
  for (int i = 0; i < Ni; ++i) {
    xs(nodes_map_inv[i]) = xy_int(i, 0);
    ys(nodes_map_inv[i]) = xy_int(i, 1);
  }
  // ADDED FOR ANALYIS
  // Output for analysis
  A_int_out = A_int, rhs_out = rhs;
}
/* SAM_LISTING_END_4 */
