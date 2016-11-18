//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <vector>
#include <complex>
#include <iostream>
#include <iomanip>
#include "timer.h"

# include <figure/figure.hpp>

//! @brief Evaluate a polynomial and its derivative using Horner scheme
//! @param[in] vector c of size $n$, coefficients of the polynomial p
//! @param[in] double x, where the polynomial has to be evaluated
//! @param[out] pair containing p(x),p'(x)
template <typename CoeffVec>
std::pair<double, double> evaldp (const CoeffVec& c, const double x) {
  std::pair<double, double> p;
  double px, dpx;
  // TODO: compute pair p using Horner scheme
  return p;
}

//! @brief Evaluate a polynomial and its derivative using a naive implementation
//! @param[in] vector c of size $n$, coefficients of the polynomial p
//! @param[in] double x, where the polynomial has to be evaluated
//! @param[out] pair containing p(x),p'(x)
template <typename CoeffVec>
std::pair<double, double> evaldp_naive(const CoeffVec& c, const double x) {
  std::pair<double, double> p;
  double px,dpx;
  // TODO: compute pair p using "naive" implementation
  return p;
}

int main() {
  std::vector<double> c {3, 1, 5, 7, 9};
  double x = .123;
    
  // Check implementations
  std::pair<double, double> p, p_naive;
  p = evaldp(c,x);
  std::cout << "Using horner scheme:\n"
            << "p(x) = " << p.first << 
            ", dp(x) = " << p.second << "\n";
    
  p_naive = evaldp_naive(c,x);
  std::cout << "Using monomial approach:\n"
            << "p(x) = " << p_naive.first 
            << ", dp(x) = " << p_naive.second << "\n";
    
  // Compare runtimes
  const unsigned repeats = 10;
   
  std::cout << std::setw(10) << "n" << std::setw(25) << "Horner scheme:" 
            << std::setw(25) << "Monomial approach:" << "\n"
            << " ================================================================\n";

  std::vector<double> e, h, m;
  for (unsigned k = 2; k <= 20; ++k) {
    Timer tm_slow, tm_fast;
    std::vector<double> c;
        
    const int n = std::pow(2, k);
    for (int i = 0; i < n; ++i) { c.push_back(i+1); }

    for (unsigned r = 0; r < repeats; ++r) {
      tm_slow.start();
      p_naive = evaldp_naive(c,x);
      tm_slow.stop();
                
      tm_fast.start();
      p = evaldp(c,x);
      tm_fast.stop();
    }
    
    std::cout << std::setw(10) << n << std::setw(25) << tm_fast.mean() 
              << std::setw(25) << tm_slow.mean() << "\n";
  }

  return 0;
}
