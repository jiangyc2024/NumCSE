//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include "autofocus.hpp"

#include <fstream>
#include <memory>

#include "pgm.hpp"
#include "conv.hpp"

// TODO: should be hidden

std::auto_ptr<Eigen::MatrixXd> M;

Eigen::MatrixXd set_focus(double f) {

    if(M.get() == nullptr) {
        std::ifstream file("image.pgm");

        PGMObject p;
        file >> p;

        M = std::auto_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd());
        p.get_data(*M);
    }

    double f0 = (1 << 1);

    double ep = std::max(std::abs(f - f0),
                         std::numeric_limits<double>::epsilon());

    unsigned int s = 16;
    Eigen::MatrixXd S(s,s);
    for(unsigned int i = 0; i < s; ++i) {
        for(unsigned int j = 0; j < s; ++j) {
            S(i,j) = 1. / (1 + (i*i + j*j)/ep);
        }
    }

    S /= S.sum();

    return conv2(*M, S);
}