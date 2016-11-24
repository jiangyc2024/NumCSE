# include <vector>
# include <Eigen/Dense>

using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
// Adaptive multilevel quadrature of a function passed in \texttt{f}.
// The vector \texttt{M} passes the positions of current quadrature nodes
template <class Function>
double adaptquad(Function& f, VectorXd& M,double rtol,double atol) {
  const std::size_t n = M.size(); // number of nodes
  VectorXd h = M.tail(n-1)-M.head(n-1), // distance of quadature nodes \Label[line]{aq:1}
           mp = 0.5*(M.head(n-1)+M.tail(n-1)); // midpoints \Label[line]{aq:2}
  // Values of integrand at nodes and midpoints \Label[line]{aq:3}
  VectorXd fx(n), fm(n - 1);
  for (unsigned i = 0; i < n; ++i) fx(i) = f(M(i));
  for (unsigned j = 0; j < n - 1; ++j) fm(j) = f(mp(j)); 
  // trapezoidal rule \eqref{eq:comptrap} \Label[line]{aq:4}
  const VectorXd trp_loc = 1./4*h.cwiseProduct(fx.head(n-1)+2*fm+fx.tail(n-1));
  // Simpson rule \eqref{eq:compsimp} \Label[line]{aq:5}
  const VectorXd simp_loc = 1./6*h.cwiseProduct(fx.head(n-1)+4*fm+fx.tail(n-1));

  // Simpson approximation for the integral value\Label[line]{aq:6}
  double I = simp_loc.sum(); 
  // local error estimate \eqref{eq:aqest}\Label[line]{aq:7}
  const VectorXd est_loc = (simp_loc-trp_loc).cwiseAbs(); 
  // estimate for quadrature error \Label[line]{aq:8}
  const double err_tot = est_loc.sum(); 

  // \com{STOP}: Termination based on \eqref{eq:aqstop}
  if ( err_tot > rtol*std::abs(I) && err_tot > atol ) { // \Label[line]{aq:term}
    // find cells where error is large
    std::vector<double> new_cells;
    for (unsigned i = 0; i < est_loc.size(); ++i) {
      // \com{MARK} by criterion \eqref{eq:mark} \& \com{REFINE} by \eqref{eq:aqref}
      if (est_loc(i) > 0.9/(n-1)*err_tot) { 
        // new quadrature point = midpoint of interval with large error
        new_cells.push_back(mp(i));
      }}

    // create new set of quadrature nodes
    // (necessary to convert std::vector to Eigen vector)
    Eigen::Map<VectorXd> tmp(new_cells.data(), new_cells.size()); 
    VectorXd new_M(M.size() + tmp.size()); 
    new_M << M, tmp; // concatenate old cells and new cells
    // nodes of a mesh are supposed to be sorted
    std::sort(new_M.data(),new_M.data()+new_M.size()); 
    I = adaptquad(f, new_M, rtol, atol); // recursion \Label[line]{aq:10}
  }
  return I;
}
/* SAM_LISTING_END_0 */
