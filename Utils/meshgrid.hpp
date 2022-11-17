# pragma once // include guard

# include <Eigen/Dense>

//! Generates a mesh, just like Matlab's meshgrid
//  Template specialization for column vectors (Eigen::VectorXd)
//  in : x, y column vectors 
//       X, Y matrices, used to save the mesh
template <typename Scalar>
void meshgrid(const Eigen::Matrix<Scalar, -1, 1>& x, 
              const Eigen::Matrix<Scalar, -1, 1>& y,
              Eigen::Matrix<Scalar, -1, -1>& X,
              Eigen::Matrix<Scalar, -1, -1>& Y) {
  const Eigen::Index nx = x.size();
  const Eigen::Index ny = y.size();
  X.resize(ny, nx);
  Y.resize(ny, nx);
  for (Eigen::Index i = 0; i < ny; ++i) {
    X.row(i) = x.transpose();
  }
  for (Eigen::Index j = 0; j < nx; ++j) {
    Y.col(j) = y;
  }
}

//! Generates a mesh, just like Matlab's meshgrid
//  Template specialization for row vectors (Eigen::RowVectorXd)
//  in : x, y row vectors 
//       X, Y matrices, used to save the mesh
template <typename Scalar>
void meshgrid(const Eigen::Matrix<Scalar, 1, -1>& x, 
              const Eigen::Matrix<Scalar, 1, -1>& y,
              Eigen::Matrix<Scalar, -1, -1>& X,
              Eigen::Matrix<Scalar, -1, -1>& Y) {
  Eigen::Matrix<Scalar, -1, 1> xt = x.transpose();
  Eigen::Matrix<Scalar, -1, 1> yt = y.transpose();
  meshgrid(xt, yt, X, Y);
}
