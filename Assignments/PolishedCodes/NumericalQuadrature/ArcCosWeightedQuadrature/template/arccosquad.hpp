
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// The function gaussquad() for pre-computed (and pre-compiled)
// Gauss-Legendre quadrature rules up to 256 nodes.
#include "gaussquad.hpp"

using namespace Eigen;

void testConvGaussQuad() {
  // Table header
  // Customize the last column or add new columns as needed.
  std::cout << std::setw(5) << "n" << std::setw(20) << " I_n " << std::setw(30)
            << " Error " << std::setw(10) << (" Behaviour\n " + std::string(65, '-'))
            << std::endl;

  // Value obtained by overkill quadrature as reference value
  constexpr double ref_val = 1.7576137811123187;

  // TO DO: Print a table that allows you to predict the asymptotic
  // behaviour of Gauss-Legendre numerical quadrature when approximating I(f).
  // START
  
  // END
}


template <typename FUNCTION>
double arccosWeightedQuad(FUNCTION &&f, unsigned int n) {
  double I = 0.0; // For accumulating quadrature result
  // TO DO: Approximate I(f) with exponential convergence in n.
  // START
  
  // END
  return I;
}

void testConvTrfGaussQuad() {
  // Table header
  // Customize the last column or add new columns as needed.
  std::cout << std::setw(5) << "n" << std::setw(20) << " I_n " << std::setw(30)
            << " Error " << std::setw(10)
            << (" Behaviour\n " + std::string(65, '-')) << std::endl;
            
  // Value obtained by overkill quadrature as reference value
  const double ref_val = 1.7576137811123187;
  
  // TO DO: Print a table that allows you to predict the asymptotic
  // behaviour of arccosWeightedQuad().
  // START
  
  // END
}

