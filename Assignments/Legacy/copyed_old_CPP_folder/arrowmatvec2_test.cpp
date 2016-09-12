#include <iostream>
#include <Eigen/Dense>
#include "arrowmatvec2.hpp"

using namespace Eigen;

// test the implementation of arrowmatvec2 with random vectors of size 5
int main(void)
{
    VectorXd a=VectorXd::Random(5);
    VectorXd d=VectorXd::Random(5);
    VectorXd x=VectorXd::Random(5);
    VectorXd Ax(5);
    
    Atimesx(d,a,x,Ax);
    VectorXd AAx(5);
    Atimesx(d,a,Ax,AAx);
    std::cout << "A*A*x = " << AAx << std::endl;
}
