/* Demonstration code for course Numerical Methods  for CSE, ETH Zurich
   Use of an object to record internal data
   @author Ralf Hiptmair
   @date September 2021
*/

#include <functional>

/* SAM_LISTING_BEGIN_9 */
template <typename RECORDER = std::function<void(int, int)>>
unsigned int myloopfunction(
    unsigned int n, unsigned int val = 1,
    RECORDER &&rec = [](int, int) -> void {}) {
  for (unsigned int i = 0; i < n; ++i) {
    rec(i, val); // Removed by the compiler for the default argument
    if (val % 2 == 0) {
      val /= 2;
    } else {
      val *= 3;
      val++;
    }
  }
  rec(n, val);
  return val;
}
/* SAM_LISTING_END_9 */
