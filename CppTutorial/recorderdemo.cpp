/* Demonstraction code for course Numerical Methods  for CS(E), ETH Zurich
   Use of a recorder object
   @author Ralf Hiptmair
   @date August 2020
*/

#include "recorderdemo.h"

#include <iostream>
#include <utility>
#include <vector>

int main(int /*argc*/, char** /*argv*/) {
  // Run without recorder
/* SAM_LISTING_BEGIN_1 */
  std::cout << "myloopfunction(10, 1) = " << myloopfunction(10, 1) << std::endl;
  // Run with recorder
  std::vector<std::pair<int, int>> store{};
  std::cout << "myloopfunction(10, 1) = "
            << myloopfunction(10, 1,
                              [&store](int n, int val) -> void {
                                store.emplace_back(n, val);
                              })
            << std::endl;
  std::cout << "History:" << std::endl;
  for (auto& i : store) {
    std::cout << i.first << " -> " << i.second << std::endl;
  }
/* SAM_LISTING_END_1 */
  return 0;
}
