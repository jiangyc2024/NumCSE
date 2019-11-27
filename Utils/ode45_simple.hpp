/* SAM_LISTING_BEGIN_0 */
template <class StateType,
          class RhsType = std::function<StateType(const StateType &)>>
class ode45 {
public:
  // Idle constructor
  ode45(const RhsType & rhs) : f(rhs) {}
  // Main timestepping routine
  template<class NormFunc = decltype(_norm<StateType>)>
  std::vector< std::pair<StateType, double>>
  solve(const StateType & y0,double T,const NormFunc & norm = _norm<StateType>);
  // Print statistics and options of this class instance.
  void print();
  struct Options {
    bool save_init = true;   // Set true if you want to save the initial data
    bool fixed_stepsize = false;  // TODO: Set true if you want a fixed step size
    unsigned int max_iterations    = 5000; 
    double min_dt = -1.; // Set the minimum step size (-1 for none)
    double max_dt = -1.; // Set the maximum step size (-1 for none)
    double initial_dt = -1.; // Set an initial step size
    double start_time = 0; // Set a starting time
    double rtol = 1e-6; // Relative tolerance for the error.
    double atol = 1e-8; // Absolute tolerance for the error.
    bool do_statistics = false; // Set to true before solving to save statistics
    bool  verbose = false;  // Print more output.
  } options;
  // Data structure for  usage statistics.
  // Available after a call of solve(), if do\_statistics is set to true.
  struct Statistics {
    unsigned int cycles = 0; // Number of loops (sum of all accepted and rejected steps)
    unsigned int steps = 0;  // Number of actual time steps performed (accepted step)
    unsigned int rejected_steps = 0; // Number of rejected steps per step
    unsigned int funcalls = 0;     // Function calls
  } statistics;
};
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
template <class StateType, class RhsType>
template <class NormFunc>
std::vector<std::pair<StateType, double>>
ode45<StateType, RhsType>::solve(const StateType &y0, double T,
                                 const NormFunc &norm) {
  const double epsilon = std::numeric_limits<double>::epsilon();
  // Setup step size default values if not provided by user
  t = options.start_time;
  // Vector for returning solution \Blue{$(t_k,\Vy_k)$}
  std::vector<std::pair<StateType, double>> snapshots;
  // The desired (initial) timestep size
  double dt = options.initial_dt;
  // Push initial data
  if (options.save_init) snapshots.push_back(std::make_pair(y0, t));
  // Temporary containers
  StateType ytemp0 = y0, ytemp1 = y0, ytemp2 = y0;
  // Pointers forswapping of temporary containers
  StateType *yprev = &ytemp0, *y4 = &ytemp1, *y5 = &ytemp2;
  // Increments \Blue{$\Vk_i$}
  std::vector<StateType> mK;  mK.resize(_s);
  // Usage statistics
  unsigned int iterations = 0; // Iterations for current step

  // \textbf{Main loop}, exit if dt too small or final time reached
  while (t < T && dt >= options.min_dt) {
    // Force hitting the endpoint of the time slot exactly
    if (t + dt > T) dt = T - t;
    // Compute the Runge-Kutta increments using the
    // coefficients provided in _mA, _vb, _vc
    mK.front() = f(*yprev);
    for (unsigned int j = 1; j < _s; ++j) {
      mK.at(j) = *yprev;
      for (unsigned int i = 0; i < j; ++i) {
        mK.at(j) += (dt * _mA(j, i)) * mK.at(i);
      }
      mK.at(j) = f(mK.at(j));
    }
    statistics.funcalls += _s;
    // Compute the 4th and the 5th order approximations
    *y4 = *yprev;
    *y5 = *yprev;
    for (unsigned int i = 0; i < _s; ++i) {
      *y4 += (dt * _vb4(i)) * mK.at(i);
      *y5 += (dt * _vb5(i)) * mK.at(i);
    }

    double tau = 2., delta = 1.;
    // Calculate the absolute local truncation error and the acceptable  error
    if (!options.fixed_stepsize) { // if (!fixed_stepsize)
      delta = norm(*y5 - *y4); // estimated 1-step error \Blue{$\mathtt{EST}_k$}
      tau = std::max(options.rtol * norm(*yprev), options.atol);
    }
    // Check if step is \com{accepted}, if so, advance
    if (delta <= tau) {
      t += dt;
      snapshots.push_back(std::make_pair(*y5, t));
      std::swap(y5, yprev);
      ++statistics.steps;
      statistics.rejected_steps += iterations;
      iterations = 0;
    }
    // Update the step size for the next integration step
    if (!options.fixed_stepsize) {
      if (delta <= std::numeric_limits<double>::epsilon()) {
        dt *= 2;
      } else {
	// Apply formula \eqref{eq:ssc}
        dt *= 0.8 * std::pow(tau / delta, _pow);
      }
      dt = std::min(options.max_dt, dt);
    } 
    ++iterations;
    ++statistics.cycles;
  } // end main loop

  return snapshots;
}
/* SAM_LISTING_END_1 */
