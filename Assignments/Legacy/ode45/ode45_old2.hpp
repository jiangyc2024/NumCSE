// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Emulates http://ch.mathworks.com/help/matlab/ref/ode45.html
// Ported from ode45.py: https://github.com/rngantner/
// See also: http://sourceforge.net/p/octave/odepkg/
 
#pragma once

#include<Eigen/Dense>

#include <iostream>
#include <limits>
#include <vector>
#include <algorithm>
#include <utility>

//! Print debugging info (comment to go silent)
// #define DEBUG

//! \brief Class for the configuration of ode45.
//! This class can be optionally passed to the ode45 function and override
//! the default selected configurations.
struct OdeConf {
    bool         save_init         = true;  //!< Set true if you want to save the initial data
    bool         fixed_stepsize    = false; //!< Set true if you want a fixed step size
    unsigned int max_iterations    = 5000;  //!< Set the maximum number of rejected iterations
    double       min_dt            = -1.;   //!< Set the minimum step size
    double       max_dt            = -1.;   //!< Set the maximum step size
    double       initial_dt        = -1.;   //!< Set an initial step size
    double       start_time        = 0;     //!< Set a starting time
};

//! \brief Class containing statistics for the ode45 solver.
//! This class is used to handle statistics for the ode45 solver.
//! Class is just a dummy/empty class for now
struct OdeStats {
    // TODO: class namaging ode45 statistics
};

//! This namespace contains definitions allowing ode45 to work with fundamental types and Eigen::Vectors.
//! Those routines are not intended for external usage.
//! Contains some template magic allowing usage of ode45 with fundamental types as StateSpace.
//! Contains magic and template trickery using SFINAE
//! DANGER: DO NOT look at this if you do not know what you are doing
namespace _private_ode45 {

//! \brief Return dimension of the vector, provided it has a Index member type and a size() member funtion
template <class T,
          typename std::enable_if<!std::numeric_limits<T>::is_specialized, bool>::type = true >
inline typename T::Index _size(const T & t) {
    return t.size();
}

//! \brief Return dimension of fundamental types as vectors, which is one
template <class T,
          typename std::enable_if<std::numeric_limits<T>::is_specialized, bool>::type = false >
inline T _size(const T & t) {
    (void) t;
    return 1;
}

//! \brief return the \f$ l^\infty \f$ norm of an Eigen::Vector.
//! T must have a member method lpNorm as Eigen::Vector
template <class T,
          typename std::enable_if<!std::numeric_limits<T>::is_specialized, bool>::type = true >
inline typename T::Scalar _norm(const T & t) {
    return t.template lpNorm<Eigen::Infinity>();
}

//! \brief return the \f$ l^\infty \f$ norm of a scalar, which is the absolute value.
template <class T,
          typename std::enable_if<std::numeric_limits<T>::is_specialized, bool>::type = false >
inline T _norm(const T & t) {
    return std::abs(t);
}

//! \brief dummy identity map if conversion is from Eigen::Matrix to Eigen::Matrix
template <class C, class T,
          typename std::enable_if<!std::numeric_limits<C>::is_specialized, bool>::type = true >
// template <class C, class T>
inline C _toStateSpace(T&& t) {
    return t;
}

//! \brief converts an Eigen::Matrix to a scalar if the matrix is 1x1
template <class C, class T,
          typename std::enable_if<std::numeric_limits<C>::is_specialized, int>::type = false >
// template <class C, class T>
inline C _toStateSpace(T&& t) {
    return t(0);
}

//! \brief dummy identity map if conversion is from Eigen::Matrix to Eigen::Matrix
template <class T,
          typename std::enable_if<!std::numeric_limits<T>::is_specialized, bool>::type = true >
inline T _toVector(T&& t) {
    return t;
}

//! \brief converts an Eigen::Matrixto a scalar if the matrix is 1x1
template <class T, typename
          std::enable_if<std::numeric_limits<T>::is_specialized, bool>::type = false >
inline Eigen::Matrix<T, 1, 1> _toVector(T&& t) {
    Eigen::Matrix<T, 1, 1> ret;
    ret(0) = t;
    return ret;
}
}

//! \brief Computes a 4-5 RK integration similar to MATLAB ode45 integrator.
//! Computes the solution of an ODE with an adaptive runge kutta method of order 4/5,
//! for the IVP \f$ y' = f(y), y(, approximated using the difference between the RK4 and RK5 methods.
//! \tparam RhsType type of the r.h.s. function \f$ f \f$, providing StateType operator()(const StateType & y)
//! \tparam StateType type of the initial data \f$ y_0 \f$ and of the solution \f$ y(t) \f$.
//! \param[in] f function for the computation of r.h.s. (e.g. a lambda function).
//! \param[in] T final time for the integration.
//! \param[in] y0 initial data \f$ y_0 = y(0) \f$.
//! \param[in] rtol relative tolerance for the error.
//! \param[in] atol absolute tolerance for the error.
//! \param[in] odeconf class containing configuration parameters for the integrations.
//! \return vector of pairs \f$ (y(t), t) \f$ for the selected snapshot times
template <class RhsType, class StateType>
std::vector< std::pair<StateType, double> > ode45(const RhsType &f,
                                                  double T, const StateType & y0,
                                                  double rtol = 1e-6, double atol = 1e-8,
                                                  OdeConf odeconf = OdeConf() ) {

    ////////////////////////////////////////
    //// 20071016, reported by Luis Randez
    // The Runge-Kutta-Fehlberg 4(5) coefficients
    // Coefficients proved on 20060827
    // See p.91 in Ascher & Petzold
    ////////////////////////////////////////
    const double       _pow = 1. / 5;
    // Number of stages
    const unsigned int _s   = 6;
    // Matrix A from the Butcher scheme
    Eigen::MatrixXd _mA(_s,_s-1);
    // Quadrature weight vectors
    Eigen::VectorXd _vb4(_s), _vb5(_s), _vc(_s);
    _mA <<  0,          0,           0,           0,            0,
            1./4.,      0,           0,           0,            0,
            3./32.,     9./32.,      0,           0,            0,
            1932./2197, -7200./2197, 7296./2197,  0,            0,
            439./216,   -8,          3680./513,   -845./4104,   0,
            -8./27.,    2,           -3544./2565, 1859./4104,   -11./40.;
    // 4th and 5th order b-coefficients for same matrix A
    _vb4 << 25./216,    0,           1408./2565,  2197./4104,   -1./5,     0;
    _vb5 << 16./135,    0,           6656./12825, 28561./56430, -9./50,    2./55;
    // Non autonomous ODEs c coefficients
    _vc = _mA.rowwise().sum();
    
    // TODO: non-autonomous ODE
    
    // Setup step size default values if not provided by user
    double t = odeconf.start_time;
    unsigned int default_nsteps = 100;
    unsigned int default_minsteps = 10;
    if(odeconf.initial_dt == -1.) {
        odeconf.initial_dt = (T - t) / default_nsteps;
    }
    if(odeconf.max_dt == -1.) {
        odeconf.max_dt = (T - t) / default_minsteps;
    }
    if(odeconf.min_dt == -1.) {
        odeconf.min_dt = (T - t) * std::numeric_limits<double>::epsilon();
    }
    
    // Solution container (returned)
    std::vector< std::pair<StateType, double> > snapshots;
    
    // Read options from odeconf
    double dt = odeconf.initial_dt;
    assert( dt > 0 && "Invalid option, dt must be positive");
    // TODO: allow negative time direction
    
    // Push initial data
    if( odeconf.save_init ) {
        snapshots.push_back( std::make_pair(y0, t) );
    }
    
    // Configuration if fixed timestepping was requested
    if( odeconf.fixed_stepsize ) {
        // TODO: implement fixed timestepping
    }
    
    // Temporary containers
    auto ytemp0 = y0;
    auto ytemp1 = y0;
    auto ytemp2 = y0;
    // Pointers forswapping of temporary containers
    auto *yprev = &ytemp0;
    auto *y4 = &ytemp1, *y5 = &ytemp2;
    
    // Dimensionality of the space
    unsigned int d = _private_ode45::_size(y0);

    // Increment matrix (reserve space)
    Eigen::MatrixXd mK(d, _s);
    
    // Usage statistics
    unsigned int iterations = 0; // Iterations for current step
    unsigned int cycles = 0;     // Number of loops
    unsigned int steps = 0;      // Number of actual steps
    
    // Main loop, exit if dt too small or final time reached
    while (t < T && dt >= odeconf.min_dt) {
        // Force hitting the endpoint of the time slot exactly
        if(t + dt > T) {
            dt = T - t;
        }

        // Compute the matrix increments using the coefficients provided in _mA, _vb, _vc
        mK.leftCols<1>() = _private_ode45::_toVector( f(*yprev) );
        for(unsigned int j = 1; j < _s; ++j) {
            mK.col(j) = _private_ode45::_toVector( f(*yprev + dt * _private_ode45::_toStateSpace<StateType>(mK.leftCols(j) * _mA.row(j).head(j).transpose()) ));
        }
        if( iterations >= odeconf.max_iterations ) {
            std::cerr << "";
            abort();
        }

        // Compute the 4th and the 5th order estimations
        *y4 = *yprev + dt * _private_ode45::_toStateSpace<StateType>(mK * _vb4);
        *y5 = *yprev + dt * _private_ode45::_toStateSpace<StateType>(mK * _vb5);

        // TODO: abs control
        double tau = 2., delta = 1.;
        // Calculate the absolute local truncation error and the acceptable  error
        if( !odeconf.fixed_stepsize) { // if (!odeconf.fixed_stepsize)
            delta = _private_ode45::_norm(*y5 - *y4);
            tau = std::max(rtol * std::max(_private_ode45::_norm(*yprev), 1.), atol);
        }
        
        // Check if step is accepted, if so, advance
        if( delta <= tau ) {
            t += dt;
            snapshots.push_back( std::make_pair(*y5, t) );
            std::swap(y5, yprev);
            ++steps;
            iterations = 0;
        }
        
        // Update the step size for the next integration step
        if( !odeconf.fixed_stepsize) {
            if ( delta <= std::numeric_limits<double>::epsilon() ) {
                dt *= 2;
            } else {
#ifdef DEBUG
                std::cout << dt << " " << tau << " " << delta << " " << 0.8 * std::pow(tau / delta, _pow) << std::endl;
#endif
                dt *= 0.8 * std::pow(tau / delta, _pow);
            }
            dt = std::min(odeconf.max_dt, dt);
        } else {
            // TODO: fixed step size
        }
        
        ++iterations;
        ++cycles;
        
        // Check if maximum number of iterations have been exeded
        if( iterations >= odeconf.max_iterations ) {
            std::cerr << "Error: Solving has not been successful. The iterative integration "
                         "loop exited at time t = " << t << " before endpoint at tend = "
                         << T << " was reached. This happened because the iterative "
                         "integration loop does not find a valid solution at this time "
                         "stamp. Try to reduce the value of \"initial_dt\" and/or \"max_dt\""
                         << std::endl;
            abort();
        }
        
    }
    
    // Check if there was a premature exit
    if( t < T) {
        std::cerr << "Error: Solving has not been successful. The iterative integration loop "
                     "exited at time t = " << t << " before endpoint at tend = " << T
                     << " was reached. This may happen if the stepsize grows smaller than "
                     "defined in \"min_dt\". Try to reduce the value of \"initial_dt\" and/or \"max_dt\""
                     << std::endl;
        abort();
    }
    
    // Returns all snapshots collected
    // TODO: option to select which snapshot to save
    return snapshots;
}
                                       
