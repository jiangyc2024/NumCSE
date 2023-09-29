/***********************************************************************
 *                                                                     *
 * Demo code                                                           *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: R.H.                                                        *
 * Date: September 2023                                                *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

// Header for basic IO
#include <functional>
#include <iostream>
#include <vector>

/* SAM_LISTING_BEGIN_1 */
int main() {
  int n_calls = 0;
  std::function<int(int)> factorial  = [&factorial, &n_calls](int n) -> int {
    n_calls++;
    if (n == 0) return 1;
    return n * factorial(n - 1);
};
  std::cout << "10! = = " << factorial(10) << std::endl;
  return 0;
}
/* SAM_LISTING_END_1 */
