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
using std::cerr;
using std::endl;

/* SAM_LISTING_BEGIN_0 */
// Auxiliary function: default norm for an \eigen vector type
template <class State>
double _norm(const State &y) { return y.norm(); } 

template <class DiscEvolOp, class State, class NormFunc = decltype(_norm<State>)>
std::vector<std::pair<double, State> >
odeintssctrl(DiscEvolOp &Psilow, unsigned int p, DiscEvolOp &Psihigh,
           State& y0,double T, double h0,
           double reltol, double abstol, double hmin,
           NormFunc &norm = _norm<State>) {
  double t = 0; // initial time \Label[line]{odeintadapt:1}
  State y = y0; // initial state
  double h = h0; // timestep to start with
  std::vector<std::pair<double,State>> states; // for output
  states.push_back({t, y});

  // Main timestepping loop
  while ((states.back().first < T) && (h >= hmin)) { // \Label[line]{ssctrl:2}
    State yh = Psihigh(h, y0); // high order discrete evolution \Blue{$\widetilde{\Psibf}^h$}\Label[line]{ssctrl:3}
    State yH = Psilow(h, y0); // low order discrete evolution \Blue{${\Psibf}^h$}\Label[line]{ssctrl:4}
    double est = norm(yH-yh); // $\leftrightarrow$ \Blue{$\mathrm{EST}_k$}\Label[line]{ssctrl:5}
    double tol = max(reltol*norm(y0), abstol); // effective tolerance \Label[line]{ssctrl:6a}

    // Optimal stepsize according to \eqref{eq:ssc}
    h = h*max(0.5,min(2.,pow(tol/est,1./(p+1)))); // \Label[line]{ssctrl:6b}
    if (est < tol)  // step \Magenta{accepted} \Label[line]{ssctrl:7}\Label[line]{ssctrl:6}
      states.push_back({t = t+min(T-t,h),y0 = yh}); // store next approximate state 
  }
  if (h < hmin) {
    cerr << "Warning: Failure at t=" 
	 << states.back().first 
	 << ". Unable to meet integration tolerances without reducing the step"
	 << " size below the smallest value allowed ("<< hmin <<") at time t." << endl;
  }
  return states;
}
/* SAM_LISTING_END_0 */
