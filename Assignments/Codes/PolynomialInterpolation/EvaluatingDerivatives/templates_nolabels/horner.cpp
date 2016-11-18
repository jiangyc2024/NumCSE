//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include "timer.h"

# include <figure/figure.hpp>

/*
 * @brief Evaluate a polynomial and its derivative using Horner scheme
 * @param[in] vector c of size $n$, coefficients of the polynomial p
 * @param[in] double x, where the polynomial has to be evaluated
 * @param[out] pair containing p(x),p'(x)
 */
template <typename CoeffVec>
std::pair<double, double> evaldp (const CoeffVec& c, const double x) {
  std::pair<double, double> p;
  double px, dpx;
  int s = c.size();

    // TODO: evaluate a polynomial using Horner scheme

  return p;
}

/*
 * @brief Evaluate a polynomial and its derivative using a naive implementation
 * @param[in] vector c of size $n$, coefficients of the polynomial p
 * @param[in] double x, where the polynomial has to be evaluated
 * @param[out] pair containing p(x),p'(x)
 */
template <typename CoeffVec>
std::pair<double, double> evaldp_naive(const CoeffVec& c, const double x) {
  std::pair<double, double> p;
  double px,dpx;
  int n=c.size();
    
    // TODO: evaluate a polynomial using naive implementation

  return p;
}

int main() {
  std::vector<double> c {3, 1, 5, 7, 9};
  double x = .123;
    
    // TODO: check implementations and compare runtimes
}
