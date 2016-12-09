//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <vector>
#include <iomanip>

#include "errors.hpp"

using namespace Eigen;

template <class Function, class State>
void rk4step(const Function &odefun, double h,
             const State & y0, State & y1)
{
    // TODO: implement a single step of the classical Runge-Kutta method of order 4
}

int main() {
    // TODO: compute $f$, $T$, $y_0$, $\VA$, $\Vb$ to run "errors(f, T, y0, A, b);"
}
