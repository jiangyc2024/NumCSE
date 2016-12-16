//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include "rkintegrator.hpp"

//! \file stabrk.cpp Solve prey/predator model with RK-SSM method

Eigen::Vector2d predprey(Eigen::Vector2d y0, double T, unsigned N)
{
    double h = T / N;
    Eigen::Vector2d y_ = y0;

    // TODO: solve the predator/prey model using the 3-stage RK-SSM of the exercise sheet

    return y_;
}

int main() {
    // TODO: solve the predator/prey model using class "RKIntegrator"
}
