# ifndef TIMER_HPP
# define TIMER_HPP

# include <chrono>
# include <iostream>
# include <vector>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * USAGE: Timer t;                                           *
 *        t.start();                                         *
 *        for (bla) { ... stuff happening ...; t.lap(); }    *
 *        double min = t.min(),                              *
 *               mean = t.mean();                            *
 *                                                           *
 *        OR                                                 *
 *                                                           *
 *        Timer t;                                           *
 *        t.start();                                         *
 *        ...  stuff happening ...                           *
 *        t.stop();                                          *
 *                                                           *
 *        NOTE: stop() and lap() are equivalent!             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

class Timer {
  
  using prec = std::chrono::nanoseconds;
  using clock = std::chrono::high_resolution_clock;
  using duration_t = std::chrono::duration<double>;
  
  private:

    static const unsigned divisor = 1e9;
    clock::time_point t_start, t_end;
    duration_t t_min = {};
    std::vector<duration_t> t_laps;

  public:

    Timer();
    void start();
    void stop();
    void lap();
    void reset();

    [[nodiscard]] double duration() const;
    [[nodiscard]] double mean() const;
    [[nodiscard]] double min() const;
};

// start the timer
inline void Timer::start(){
  t_start = clock::now();
  t_end = t_start;
}

// stop the timer (equivalent to lap)
inline void Timer::stop(){
  // stop is just another lap
  lap();
}

// new lap
inline void Timer::lap(){
  const clock::time_point tmp = clock::now();

  // get laptime
  duration_t laptime;
  laptime = tmp - t_end;

  // check if this lap was faster
  if (t_min > laptime || t_laps.empty()) {
    t_min = laptime;
  }

  // save time of this lap
  t_laps.push_back(laptime);

  // save total time
  t_end = tmp;
}

// idle constructor
inline Timer::Timer() {
  #ifndef NDEBUG
  static bool runonce = true;
  if (runonce) {
    std::cerr << "Warning: Timer was build as DEBUG." << std::endl;
    runonce = false;
  }
  #endif
}

// resets all values
inline void Timer::reset(){
  t_laps = std::vector<duration_t>();
  start();
}

// returns total duration timer has been running
inline double Timer::duration() const {
  double result = 0.;
  if (!t_laps.empty()){
    // returning time in seconds! thats what the divisor is for
    auto dur = std::chrono::duration_cast<prec>(t_end - t_start);
    result = static_cast<double>(dur.count())/divisor;
  }
  else {
    std::cerr << "Before calling Timer::duration() you need to call Timer::lap() or Timer::stop()!\n";
  }
  return result;
}

// returns mean of all laps
inline double Timer::mean() const {
  double result = 0.;
  if (!t_laps.empty()){
    // save total time in std::chrono units
    auto total_time = t_laps[0];
    for (unsigned int i = 1; i < t_laps.size(); ++i) {
      total_time += t_laps[i];
    }
    // convert time to double
    auto total_dur = std::chrono::duration_cast<prec>(total_time);
    result = static_cast<double>(total_dur.count())/divisor/static_cast<double>(t_laps.size());
  }
  else {
    std::cerr << "Before calling Timer::mean() you need to call Timer::lap() or Timer::stop()!\n";
  }
  return result;
}

// returns minimum of all laps
inline double Timer::min() const {
  auto min = std::chrono::duration_cast<prec>(t_min);
  return static_cast<double>(min.count())/divisor;
}

# endif
