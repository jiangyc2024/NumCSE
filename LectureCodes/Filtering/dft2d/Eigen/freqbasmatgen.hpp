///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>

// Generation of Fourier basis matrix of size mxn
inline Eigen::MatrixXcd freqbasmatgen(int m,int n,int k,int l) {
  using Comp = std::complex<double>;
  const Comp iu = Comp(0.0,1.0);
  // The 2D Fourier basis matrices are tensor products of 1D Fourier basis vectors
  Eigen::VectorXcd u(m); const Comp om = exp((2.0*M_PI*iu)*(static_cast<double>(k)/m));
  Comp t=1.0;
  for (int j=0;j<m;j++,t*=om) {
    u(j) = t;
  }
  Eigen::RowVectorXcd v(n); const Comp on = exp((2.0*M_PI*iu)*(static_cast<double>(l)/n));
  t=1.0; 
  for(int j=0;j<n;j++,t*=on) {
    v(j) = t;
  }
  Eigen::MatrixXcd B(m,n); B = u*v;
  return B;
}

