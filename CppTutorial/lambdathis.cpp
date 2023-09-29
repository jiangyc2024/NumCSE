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
#include <algorithm>
#include <iostream>
#include <vector>

/* SAM_LISTING_BEGIN_2 */
struct X {
  explicit X(int N) : N_(N) {}
  unsigned int mod(unsigned int n) const { return N_ % n; }
  bool modmatch(const std::vector<int> nums, unsigned int n) const;
  int N_;
};

bool X::modmatch(const std::vector<int> nums, unsigned int n) const {
  auto it =
    std::find_if(nums.begin(), nums.end(), [this, n](int k) -> bool {
      return (N_ != 0) and ((k % n) == mod(n));
    });
  return (it != nums.end());
}

/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_1 */
int main() {
  X x(7);
  std::cout << (x.modmatch({2,11,6,8},3)?"mod contained":"mod not contained") << std::endl;
  std::cout << (x.modmatch({1,3,5,7},3)?"mod contained":"mod not contained") << std::endl;
}
/* SAM_LISTING_END_1 */
