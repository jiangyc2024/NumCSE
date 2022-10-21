#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>
#include <functional>

using namespace std;
using namespace std::chrono;

/**
 * @brief time an action by minimum execution time of N repetitions
 * @param a action (lambda) to time
 * @param N number of repetitions
 * @return minimum execution time in microseconds
 */
template <class Action>
size_t time(Action&& a, size_t const N) {
  size_t minimum = ~0;
  for (size_t i = 0; i < N; ++i) {
    auto start = high_resolution_clock::now();
    a();
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    minimum = std::min(minimum, (size_t)duration.count());
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
               const vector<function<void()>>& actions, InitF&& init,
               PostF&& post, const size_t N, const size_t w,
               const vector<string>& header = {}) {
  if (header.size()) {
    for (auto h : header) cout << setw(w) << h;
    cout << endl;
  }

  for (size_t i(0); i < params.size(); ++i) {
    init(params[i]);
    cout << setw(w) << params[i];
    for (auto a : actions) {
      cout << setw(w) << time(a, N);
    }
    cout << endl;
    post();
  }
}

template <class ParamVector, class InitF>
void timeTable(const ParamVector& params,
               const vector<function<void()>>& actions, InitF&& init,
               const vector<string>& header = {}) {
  timeTable(
      params, actions, init, [] {}, 5, 10, header);
}

template <class ParamVector, class InitF, class PostF>
void timeTable(const ParamVector& params,
               const vector<function<void()>>& actions, InitF&& init,
               PostF&& post, const vector<string>& header = {}) {
  timeTable(params, actions, init, post, 5, 10, header);
}

template <class ParamVector, class InitF>
void timeTable(const ParamVector& params,
               const vector<function<void()>>& actions, InitF&& init,
               const size_t N, const size_t w,
               const vector<string>& header = {}) {
  timeTable( params, actions, init, [] {}, N, w, header );
}
