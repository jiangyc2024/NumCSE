///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Example {ex:ode45stiff}
/// (C) 2016 SAM, D-MATH
/// Author(s): R. Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <iostream>
#include <math.h>

#define MATLABCOEFF true

#include "ode45.hpp"
#include "ode45_boost.hpp"

using namespace std;

int main(int argc,char **argv)
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << "[1-2]" << std::endl;
    exit(-1L);
  }
  else {
    const int select = atoi(argv[1]);
    switch(select) {
    case 0: {
    std::cerr << "Usage: " << argv[0] << "[1-2]" << std::endl;
    exit(-1L);
    }
    case 1: {
  cout << "Explicit adaptive RK-SSM for stiff scalar logistic IVP" << endl;
/* SAM_LISTING_BEGIN_0 */
  // Types to be used for a scalar ODE with state space \Blue{$\bbR$}
  using StateType = double;
  using RhsType = std::function<StateType(StateType)>;
  // Logistic differential equation \eqref{eq:logode}
  double lambda = 500.0;
  RhsType f = [lambda](StateType y) { return lambda*y*y*(1-y); };
  StateType y0 = 0.01;   // Initial value, will create a STIFF IVP
  // State space \Blue{$\bbR$}, simple modulus supplies norm
  auto normFunc = [](StateType x){ return fabs(x); };
  
  // Invoke explicit Runge-Kutta method with stepsize control
  ode45<StateType,RhsType> integrator(f);
  // Set rather loose tolerances
  integrator.options.rtol = 0.1;
  integrator.options.atol = 0.001;
  integrator.options.min_dt = 1E-18;
  std::vector<std::pair<StateType, double>> states = integrator.solve(y0,1.0,normFunc);
  // Output information accumulation during numerical integration
  integrator.options.do_statistics = true; integrator.print();
/* SAM_LISTING_END_0 */
  
  for (auto state : states)
    std::cout << "t = " << state.second << ", y = " << state.first << std::endl;

  // Alternative numerical integrator
  // Invoke explicit Runge-Kutta method with stepsize control
  cout << endl << "\t RKF5(4) integrator " << endl;
  std::vector<std::pair<double, double>> res = rkf45(f,y0,1.0,normFunc,0.01,0.1,0.001,1E-18);

  for (auto state : res)
    std::cout << "t = " << state.first << ", y = " << state.second << std::endl;
  break;
  }
    case 2: {
   cout << "Explicit adaptive RK-SSM for stiff decay IVP" << endl;
/* SAM_LISTING_BEGIN_1 */
  // Types to be used for a scalar ODE with state space \Blue{$\bbR$}
  using StateType = double;
  using RhsType = std::function<StateType(StateType)>;
  // Decay differential equation \eqref{eq:logode}
  double lambda = 80.0;
  RhsType f = [lambda](StateType y) { return -lambda*y; };
  StateType y0 = 1.0;   // Initial value, will create a STIFF IVP
  // State space \Blue{$\bbR$}, simple modulus supplies norm
  auto normFunc = [](StateType x){ return fabs(x); };
  
  // Invoke explicit Runge-Kutta method with stepsize control
  ode45<StateType,RhsType> integrator(f);
  // Set rather loose tolerances
  integrator.options.rtol = 0.1;
  integrator.options.atol = 0.001;
  integrator.options.min_dt = 1E-18;
  std::vector<std::pair<StateType, double>> states = integrator.solve(y0,1.0,normFunc);
  // Output information accumulation during numerical integration
  integrator.options.do_statistics = true; integrator.print();
/* SAM_LISTING_END_1 */
  
  for (auto state : states)
    std::cout << "t = " << state.second << ", y = " << state.first << std::endl;

  // Alternative numerical integrator
  // Invoke explicit Runge-Kutta method with stepsize control
  cout << endl << "\t RKF5(4) integrator " << endl;
  std::vector<std::pair<double, double>> res = rkf45(f,y0,1.0,normFunc,0.01,0.1,0.001,1E-18);

  for (auto state : res)
    std::cout << "t = " << state.first << ", y = " << state.second << std::endl;
  break;
  }}}
}
