#include <Eigen/Dense>
#include <iostream>
#include "arrowmatvec.hpp"

using namespace Eigen;

// We test the function arrowmatvec with 5 dimensional random vectors.
int main(void){
//     srand((unsigned int) time(0));
    VectorXd a=VectorXd::Random(5);
    VectorXd d=VectorXd::Random(5);
    VectorXd x=VectorXd::Random(5);
    VectorXd y;
    
    arrowmatvec(d,a,x,y);
    std::cout << "A*A*x = " << y << std::endl;
}
