#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

/**
 * @brief time an action by minimum execution time of N repetitions
 * @param a action (lambda) to time
 * @param N number of repetitions
 * @return minimum execution time in microseconds
 */
template <class Action>
size_t time(Action&& a, size_t const N) {
  size_t minimum = std::numeric_limits<size_t>::max();
  for (size_t i = 0; i < N; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
    a();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    minimum = std::min(minimum, static_cast<size_t>(duration.count()));
  }

  return minimum;
}

/**
 * @brief print a table of all actions timed over all parameters
 * @param params vector of parameters
 * @param actions vector of actions
 * @param init initialization function depending on parameter, called once for
 * each parameter
 * @param post post processing/cleanup function, called once after each timed
 * action
 * @param N number of repetitions
 * @param w column width
 * @param header (optional) table column names
 */
template <class ParamVector, class InitF, class PostF>
void timeTable(const ParamVector& params,
               const std::vector<std::function<void()>>& actions, InitF&& init,
               PostF&& post, const size_t N, const int w,
               const std::vector<std::string>& header = {}) {
  if (!header.empty()) {
    for (auto const & h : header) {
      std::cout << std::setw(w) << h;
    }
    std::cout << std::endl;
  }

  for (typename ParamVector::Index i(0); i < params.size(); ++i) {
    init(params[i]);
    std::cout << std::setw(w) << params[i];
    for (auto const & a : actions) {
      std::cout << std::setw(w) << time(a, N);
    }
    std::cout << std::endl;
    post();
  }
}

template <class ParamVector, class InitF>
void timeTable(const ParamVector& params,
               const std::vector<std::function<void()>>& actions, InitF&& init,
               const std::vector<std::string>& header = {}) {
  timeTable(
      params, actions, init, [] {}, 5, 10, header);
}

template <class ParamVector, class InitF, class PostF>
void timeTable(const ParamVector& params,
               const std::vector<std::function<void()>>& actions, InitF&& init,
               PostF&& post, const std::vector<std::string>& header = {}) {
  timeTable(params, actions, init, post, 5, 10, header);
}

template <class ParamVector, class InitF>
void timeTable(const ParamVector& params,
               const std::vector<std::function<void()>>& actions, InitF&& init,
               const size_t N, const size_t w,
               const std::vector<std::string>& header = {}) {
  timeTable( params, actions, init, [] {}, N, w, header );
}
