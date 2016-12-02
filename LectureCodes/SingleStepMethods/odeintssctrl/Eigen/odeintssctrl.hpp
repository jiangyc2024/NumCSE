///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <math.h>

using std::max;
using std::min;

/* SAM_LISTING_BEGIN_0 */
template <class Func, class Func2, class NormFunc, class State>
std::vector<std::pair<double, State>> odeintssctrl(
    Func &Psilow, unsigned int p, Func2 &Psihigh,NormFunc &norm, State& y0,
    double T, double h0, double reltol, double abstol, double hmin) {
  double t = 0; // \Label[line]{ssctrl:1}
  State y = y0;
  double h = h0;
  
  std::vector<std::pair<double,State>> states;
  states.push_back({t, y});
  
  while ((states.back().first < T) && (h >= hmin)) // \Label[line]{ssctrl:2}
    {
      State yh = Psihigh(h, y0); // high order discrete evolution \Blue{$\widetilde{\Psibf}^h$}\Label[line]{ssctrl:3}
      State yH = Psilow(h, y0); // low order discrete evolution \Blue{${\Psibf}^h$}\Label[line]{ssctrl:4}
      double est = norm(yH-yh); // $\leftrightarrow$ \Blue{$\mathrm{EST}_k$}\Label[line]{ssctrl:5}
      double tol = std::max(reltol*norm(y0), abstol); // \Label[line]{ssctrl:6a}
      
      h = h*max(0.5,min(2.,pow(tol/est,1./(p+1)))); // Optimal stepsize according to \eqref{eq:ssc}\Label[line]{ssctrl:6b}
      if (est < tol)  { // \Label[line]{ssctrl:6}
	//step \Magenta{accepted} \Label[line]{ssctrl:7}
	y0 = yh; 
	states.push_back({states.back().first+min(T-states.back().first, h), y0});
      }
    }
  if (h < hmin) {
      std::cerr << "Warning: Failure at t=" 
		<< states.back().first 
		<< ". Unable to meet integration tolerances without reducing the step size below the smallest value allowed ("<< hmin <<") at time t." << std::endl;
  }
  return states;
}
/* SAM_LISTING_END_0 */
