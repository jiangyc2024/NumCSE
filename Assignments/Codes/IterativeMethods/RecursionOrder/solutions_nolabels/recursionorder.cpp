//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <iostream>
#include <cmath>

int main(int argc, char**argv) {
    int n = 10;

    std::vector<double> e(n), loge(n);

    e[0] = 1;   loge[0] = std::log(e[0]);
    e[1] = 0.8; loge[1] = std::log(e[1]);
    for (int k=1; k < n; ++k) {
        e[k+1] = e[k]*std::sqrt(e[k-1]);
        loge[k+1] = std::log(e[k+1]);

        std::cout << (loge[k+1] - loge[k] ) / (loge[k] - loge[k-1])
                  << std::endl;
    }
}
