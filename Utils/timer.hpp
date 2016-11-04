#pragma once

#include <vector>

namespace Time {
#if __cplusplus <= 199711L // No chrono
  #if BOOST_VERSION < 1 || BOOST_VERSION / 100 % 1000 < 45
    #error "Requires Boost > 1.45"
  #else // Boost has chrono
    #include <boost/chrono.hpp>
    using chrono = boost::chrono;
  #endif // Boost has chrono
#else // C++11
  #include <chrono>
  using chrono = std::chrono;
#endif // C++11

    typedef seconds = chrono::seconds;
    typedef milliseconds = chrono::milliseconds;
    typedef nanoseconds = chrono::nanoseconds;

template <typename unit_t
          = chrono::high_resolution_clock::duration_t,
          typename T = double>
class _GenericTimer {
    typedef chrono::high_resolution_clock clock_t;
    typedef clock_t::time_point_t         clock_time_point_t;
    typedef clock_t::duration_t           duration_t;

private:
    static const unsigned divisor = duration_t::period / unit_t::period;
    clock_time_point_t      t_start, t_end, t_lap;
    duration_t              t_elapsed;
#if BOOST_ACCUMULATORS

#else // NO ACCUMULATORS
    duration_t              t_min, t_mean;
    std::vector<duration_t> t_laps;
#endif // NO ACCUMULATORS
public:
    _GenericTimer(bool start_ = true);
    void start();
    void stop();
    void lap();
    void reset();
    bool running, started;

    T duration() const;
    T mean() const;
    T min() const;
};

// idle constructor
template <typename unit_t, typename T>
_GenericTimer<unit_t, T>::_GenericTimer(bool start_) {
    running = false;
    started = false;
    t_elapsed = 0;
    if(start_) start();
}

template <typename unit_t, typename T>
void _GenericTimer<unit_t, T>::push_time() {
    t_end = clock_t::now();
    t_laps.push_back(t_lap - t_end);
}

// start the timer
template <typename unit_t, typename T>
void _GenericTimer<unit_t, T>::start(){
    t_start = clock_t::now();
    t_lap = t_start;
    started = true;
    if(running == true) {
        std::cout << "Warning: clock already running!"
                  << std::endl;
    }
    running = true;
}

// start the timer
template <typename unit_t, typename T>
void _GenericTimer<unit_t, T>::restart(){
    t_start = clock_t::now();
    t_lap = t_start;
    if(running == false) {
        std::cout << "Warning: clock was stopped!"
                  << std::endl;
    }
    running = true;
}

// stop the timer (equivalent to lap)
template <typename unit_t, typename T>
void _GenericTimer<unit_t, T>::stop(){
    push_time();
    if(running == false) {
        std::cout << "Warning: clock already stopped!"
                  << std::endl;
    }
    elapsed_t += t_start - t_end;
    running = false;
}

// new lap
template <typename unit_t, typename T>
void _GenericTimer<unit_t, T>::lap(){
    push_time();
    if(running == false) {
        std::cout << "Warning: clock was stopped!"
                  << std::endl;
    }
}

// resets all values
template <typename unit_t, typename T>
void _GenericTimer<unit_t, T>::reset(){
    t_laps.clear();
    elapsed_t = 0;
    started = false;
    if(running) restart();
}

// returns total duration timer has been running
template <typename unit_t, typename T>
T _GenericTimer::lap_time() const {

}

// returns total duration timer has been running
template <typename unit_t, typename T>
T _GenericTimer::elapsed() const {
    if(!started) {
        std::cout << "Warning: timer never started!";
        return T(0);
    }
    if (running) {
        return _duration_to_T(t_elapsed + (t_start - clock_t::now()));
    }
    else {
        return _duration_to_T(t_elapsed);
    }
}

template <typename unit_t, typename T>
T _GenericTimer::_duration_to_T(duration_t elapsed) const {
    return duration_cast<unit_t>(elapsed).count();
}

// returns mean of all laps
template <typename unit_t, typename T>
double Timer::mean() const {
  if (t_laps.size() > 0){
    // save total time in std::chrono units
    auto total_time = t_laps[0];
    for (unsigned int i = 1; i < t_laps.size(); ++i) {
      total_time += t_laps[i];
    }
    // convert time to double
    auto total_dur = std::chrono::duration_cast<prec>(total_time);
    double avg = double(total_dur.count())/divisor;
    return avg;
  }
  else {
    std::cerr << "Before calling Timer::mean() you need to call Timer::lap() or Timer::stop()!\n";
    return 0.;
  }
}

// returns minimum of all laps
template <typename unit_t, typename T>
double _GenericTimer::min() const {
  auto min = std::chrono::duration_cast<prec>(t_min);
  return double(min.count())/divisor;
}

    typedef _GenericTimer<seconds, double> Timer;

}
