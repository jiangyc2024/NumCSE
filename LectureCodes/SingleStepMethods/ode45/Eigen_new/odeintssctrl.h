///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <vector>

/* SAM_LISTING_BEGIN_0 */
// Auxiliary function: default norm for an \eigen vector type
template <class State>
double defaultnorm(const State &y) { return y.norm(); }
// Auxilary struct to hold user options
struct Odeintssctrl_options {
  double T;       // terminal time
  double h0;      // initial time step
  double reltol;  // norm-relative error tolerance
  double abstol;  // absolute error tolerance
  double hmin;    // smallest allowed time step
} __attribute__((aligned(64)));
// Adaptive single-step integrator
template <class DiscEvolOp, class State,
          class NormFunc = decltype(defaultnorm<State>)>
std::vector<std::pair<double, State>> odeintssctrl(
    DiscEvolOp &&Psilow, unsigned int p, DiscEvolOp &&Psihigh, const State &y0,
    const Odeintssctrl_options &opt,
    NormFunc &norm = defaultnorm<State>) {
  double t = 0;           // initial time $\cob{t_0=0}$\Label[line]{odeintadapt:1}
  State y = y0;           // current state, initialized here
  double h = opt.h0;  // timestep to start with
  // Array for returning pairs of times/states $\cob{\left(t_k,\Vy_k\right)_k}$
  std::vector<std::pair<double, State>> states;  
  states.push_back({t, y}); // initial time and state
  // Main timestepping loop
  while ((states.back().first < opt.T) && (h >= opt.hmin)) {  // \Label[line]{ssctrl:2}
    State yh = Psihigh(h, y);  // high-order discrete evolution \Blue{$\widetilde{\Psibf}^h$}\Label[line]{ssctrl:3}
    State yH = Psilow(h, y);  // low-order discrete evolution \Blue{${\Psibf}^h$}\Label[line]{ssctrl:4}
    const double est = norm(yH - yh);  // $\leftrightarrow$ \Blue{$\mathrm{EST}_k$}\Label[line]{ssctrl:5}
    const double tol = std::max(opt.reltol * norm(y), opt.abstol);  // effective tolerance \Label[line]{ssctrl:6a}
    // Optimal stepsize according to \eqref{eq:ssc}
    if (est < tol) {  // step \Magenta{accepted} \Label[line]{ssctrl:7}\Label[line]{ssctrl:6}
      // store next approximate state
      states.push_back({t = t + std::min(opt.T - t, h), y = yh});  
    }
    // New timestep size according to \eqref{eq:ssc}
    h *= std::max(0.5, std::min(2., 0.9 * std::pow(tol / est, 1. / (p + 1)))); //\Label[line]{ssctrl:6b}
    if (h < opt.hmin) {
      std::cerr << "Warning: Failure at t=" << states.back().first
          << ". Unable to meet integration tolerances without reducing the step"
          << " size below the smallest value allowed (" << opt.hmin
          << ") at time t = " << t << "." << std::endl;
    }
  } // end main loop
  return states;
}
/* SAM_LISTING_END_0 */
