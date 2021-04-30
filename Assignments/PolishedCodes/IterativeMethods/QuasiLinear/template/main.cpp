///
/// Minimal runner for (9-10)
///

#include <cmath>
#include "quasilin.hpp"

int main(void) {
    
    constexpr unsigned int n = 5;
    Eigen::VectorXd x = Eigen::VectorXd::Constant(n, 1 / std::sqrt(n));
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    
    // test implementation of fixed_point_step and newton_step()
    Eigen::VectorXd x_fp = fixed_point_step(x, b);
    std::cout << "Single fixed-point step: " << std::endl 
              << x_fp << std::endl << std::endl;
    
    
    Eigen::VectorXd x_nwt = newton_step(x, b);
    std::cout << "Single Newton step:" << std::endl
              << x_nwt << std::endl << std::endl;
              
    // test iterative methods         
    int testN;
    constexpr double atol = 1e-13;
    constexpr double rtol = 1e-11;
    
    std::cout << "Enter 1 to test the fixed point \nEnter 2 to test the Newton method \nEnter 0 to escape\n"; 
    std::cin >> testN;
    std::cout.precision(10);
    switch (testN) {
        case 1: x_fp = solveQuasiLinSystem(rtol, atol, b);
                std::cout << "Fixed point method converges to: " << std::endl
                          << x_fp << std::endl << std::endl;
                break;
                          
        case 2: x_nwt = solveQLSystem_Newton(rtol, atol, b);
                std::cout << "Newton method converges to:" << std::endl
                          << x_nwt << std::endl << std::endl; 
                break;
        default: break;
    }
}
