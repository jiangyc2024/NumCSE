#include <iostream>

#include <Eigen/Dense>

#include "timer.h"
#include "conv.hpp"

using namespace Eigen;

int main(int argc, char **argv) {

    if(argc <= 1) {
        MatrixXd A(5,5);
        MatrixXd B(3,3);

        A << 1,2,3,4,5,
                4,5,6,7,8,
                7,8,9,10,11,
                1,1,1,1,1,
                2,2,2,2,2;
        B << 0,1,0,
                1,2,1,
                0,1,1;

        std::cout << conv2(A, B);
    } else {
        unsigned int n = std::stoi(argv[1]);
        MatrixXd A = MatrixXd::Random(n, n);
        unsigned int m = n;
        if(argc > 2) m = std::stoi(argv[2]);
        MatrixXd B = MatrixXd::Random(m, m);

        Timer t;
        t.start();
        MatrixXd C = conv2(A,B);
        t.stop();
        std::cout << "Elapsed:" << t.duration() << "s" << std::endl;
        std::cout << "Norm:" << C.norm() << std::endl;
    }
}
