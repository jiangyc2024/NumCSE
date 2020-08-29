#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;
using namespace std::chrono;

/**
 * @brief time an action by minimum of N repetitions
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
 * @param init initialization function depending on parameter
 * @param N number of repetitions
 * @param header (optional) table column names
 */
template <class Param>
void timeTable( const vector<Param>& params, const vector<function<void()>>& actions, void (&init)(Param), const size_t N, const size_t w, const vector<string> & header = {}) {
  
  if( header.size( )) {
    for( auto h : header ) cout << setw(w) << h;
    cout << endl;
  }

  for (auto p : params) {
    cout << setw(w) << p;
    init(p);
    for (auto a : actions) cout << setw(w) << time(a, N);
    cout << endl;
  }
}

template <class Param>
void timeTable( const vector<Param>& params, const vector<function<void()>>& actions, void (&init)(Param), const vector<string> & header = {} ) {
  timeTable( params, actions, init, 5, 10, header );
}

