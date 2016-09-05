template<class VecType, class MatType>
VecType invpowit(const Eigen::MatrixBase<MatType> &A,double tol)
{
  using index_t = typename MatType::Index;
  using scalar_t = typename VecType::Scalar;
  // Make sure that the function is called with a square matrix
  const index_t n = A.cols();
  const index_t m = A.rows();
  eigen_assert(n == m);
  // Request LU-decomposition
  auto A_lu_dec = A.lu();
  // Initial guess for inverse power iteration
  VecType xo = VecType::Zero(n);
  VecType xn = VecType::Random(n); 
  // Normalize vector
  xn /= xn.norm();
  // Terminate if relative (normwise) change below threshold
  while ((xo-xn).norm() > xn.norm()*tol) {
    xo = xn;
    xn = A_lu_dec.solve(xo);
    xn /= xn.norm();
  }
  return(xn);
}
