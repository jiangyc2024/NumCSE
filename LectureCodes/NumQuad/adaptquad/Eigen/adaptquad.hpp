# include <Eigen/Dense>

using Eigen::VectorXd;

template <class Function>
double adaptquad(Function& f, VectorXd& M, const double& rtol, const double& atol) {
  const unsigned n = M.size(); // number of nodes
  VectorXd h = M.tail(n - 1) - M.head(n - 1), // distance of quadature nodes \label{aq:1}
           mp = 0.5*(M.head(n - 1) + M.tail(n - 1)); // midpoints of intervals \label{aq:2}
  
  // values at nodes and midpoints \label{aq:3}
  VectorXd fx(n), fm(n - 1);
  for (unsigned i = 0; i < n; ++i) {
    fx(i) = f(M(i));
  }
  for (unsigned j = 0; j < n - 1; ++j) {
    fm(j) = f(mp(j));
  }

  // trapezoidal rule \eqref{eq:comptrap}\label{aq:4}
  const VectorXd trp_loc = 1./4*h.cwiseProduct(fx.head(n - 1) + 2*fm + fx.tail(n - 1));
  // Simpson rule \eqref{eq:compsimp}\label{aq:5}
  const VectorXd simp_loc = 1./6*h.cwiseProduct(fx.head(n - 1) + 4*fm + fx.tail(n - 1));

  double I = simp_loc.sum(); // Simpson approximation for the integral value \label{aq:6}
  const VectorXd est_loc = (simp_loc - trp_loc).cwiseAbs(); // local error estimate \eqref{eq:aqest}\label{aq:7}
  const double err_tot = est_loc.sum(); // estimate for quadrature error \label{aq:8}

  // Termination based on \eqref{eq:aqstop}
  if ( err_tot > rtol*std::abs(I) && err_tot > atol ) { // \label{aq:term}
    // find cells where error is large
    std::vector<double> new_cells;
    for (unsigned i = 0; i < est_loc.size(); ++i) {
      // if local error > 0.9/(no of intervals)*(total error) add new cell
      if (est_loc(i) > 0.9/(n - 1)*err_tot) { 
        // new quadrature point = midpoint of interval with large error
        new_cells.push_back(mp(i));
      }
    }

    // create new mesh:
    Eigen::Map<VectorXd> tmp(new_cells.data(), new_cells.size()); // convert std::vector to Eigen vector
    VectorXd new_M(M.size() + tmp.size()); 
    new_M << M, tmp; // concatenate old cells and new cells 
    std::sort(new_M.data(), new_M.data() + new_M.size()); // sort new mesh

    I = adaptquad(f, new_M, rtol, atol); // recurse \label{aq:10}
  }
  return I;
}
