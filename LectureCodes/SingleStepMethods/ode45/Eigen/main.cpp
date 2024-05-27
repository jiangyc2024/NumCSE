///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include "ode45.hpp"
#include <iostream> 

int main()
{
  // Logistic differential equation \eqref{eq:logode}
  auto f = [](double x){ return 5*x*(1-x); };
  // State space \Blue{$\bbR$}, simple modulus supplies norm
  auto normFunc = [](double x){ return fabs(x); };
  
  const double y0 = 0.2; // initial value
  // Invoke explicit Runge-Kutta-Fehlberg method with stepsize control
  const std::vector<std::pair<double, double>> states = rkf45(f, y0, 1, normFunc);
  
  for (auto state : states) {
    std::cout << "t = " << state.first << ", y = " << state.second << std::endl;
  }
}
