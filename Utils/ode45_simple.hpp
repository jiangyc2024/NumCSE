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

