# include <Eigen/Dense>

using Eigen::VectorXd;

void divdiff(const VectorXd& t, const VectorXd& y, VectorXd& c) {
  // IN:  t = node set (mutually different)
  //      y = nodal values 
  // OUT: c = coefficients of polynomial in Newton basis
  
  c = y;
  const unsigned n = y.size() - 1;
  for (unsigned l = 0; l < n; ++l) {
    for (unsigned j = l; j < n; ++j) {
      c(j+1) = ( c(j+1) - c(l) )/( t(j+1) - t(l) );
    }
  }
}
