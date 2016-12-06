//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

//! \brief One step of autonomous IVP y'' + y' + y = 0, [y(0), y'(0)] = z0 using SDIRK method
//! Use SDIRK method for first order ode z' = f(z). Steps of size h.
//! \tparam StateType type of solution space y and initial data y0
//! \param[in] z0 initial data z(0)
//! \param[in] h size of the step
//! \param[in] gamma parameter
//! \return next step z1
template <class StateType>
StateType sdirkStep(const StateType & z0, double h, double gamma) {
    // TODO: implement one step of SDIRK method
    return z0; // Placeholder
}

//! \brief Solve autonomous IVP y'' + y' + y = 0, [y(0), y'(0)] = z0 using SDIRK method
//! Use SDIRK method for first order ode z' = f(z), with N equidistant steps
//! \tparam StateType type of solution space z = [y,y']! and initial data z0 = [y(0), y'(0)]
//! \param[in] z0 initial data z(0)
//! \param[in] N number of equidistant steps
//! \param[in] T final time of simulation
//! \param[in] gamma parameter
//! \return vector containing each step of z_k (y and y')
template <class StateType>
std::vector<StateType> sdirkSolve(const StateType & z0, unsigned int N, double T, double gamma) {
    // Solution vector
    std::vector<StateType> res(N+1);
    // Equidistant step size
    const double h = T/N;
    
    // Push initial data
    res.at(0) = z0;
    
    // Main loop
    for(unsigned int i = 1; i <= N; ++i) {
        res.at(i) = sdirkStep(res.at(i-1), h, gamma);
    }
    
    return res;
}

int main() {
    // Initial data z0 = [y(0), y'(0)]
    Eigen::Vector2d z0;
    z0 << 1,0;
    // Final time
    const double T = 10;
    // Parameter 
    const double gamma = (3.+std::sqrt(3.)) / 6.;
    // Mesh sizes
    std::vector<int> N = {20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240};
    // Exact solution (only y(t)) given z0 = [y(0), y'(0)] and t
    auto yex = [&z0] (double t) {
        return 1./3.*std::exp(-t/2.) * ( 3.*z0(0) * std::cos( std::sqrt(3.)*t/2. ) + 
        std::sqrt(3.)*z0(0) * std::sin( std::sqrt(3.)*t/2. ) + 
        2.*std::sqrt(3.)*z0(1) * std::sin( std::sqrt(3.)*t/2. ) );
    };
    
    // Store old error for rate computation
    double errold = 0;
    std::cout << std::setw(15) << "n" << std::setw(15) << "maxerr" << std::setw(15) << "rate" << std::endl;
    // Loop over all meshes
    for(unsigned int i = 0; i < N.size(); ++i) {
        int n = N.at(i);
        // Get solution
        auto sol =  sdirkSolve(z0, n, T, gamma);
        // Compute error
        double err = std::abs(sol.back()(0) - yex(T));
        
        // I/O
        std::cout << std::setw(15) << n << std::setw(15) << err;
        if(i > 0) std::cout << std::setw(15) << std::log2(errold /err);
        std::cout << std::endl;
        
        // Store old error
        errold = err;
    }
}
