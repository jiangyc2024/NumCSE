/***********************************************************************
 *                                                                     *
 * Demo code                                                           *    
 * (Prof. Dr. R. Hiptmair)                                             * 
 * Author: R.H.                                                        *
 * Date: May 2018                                                      * 
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

// Header for basic IO
#include <iostream>
// Provides random acccess container class
#include <vector>
// Provides algorithms operating on generic containers
#include <algorithm>

using namespace std;

/* SAM_LISTING_BEGIN_1 */
int main() {
  // initialize a vector from an initializer list
  std::vector<double> v({1.2,2.3,3.4,4.5,5.6,6.7,7.8});
  // A vector of the same length
  std::vector<double> w(v.size());
  // Do cumulative summation of v and store result in w
  double sum = 0;
  std::transform(v.begin(),v.end(),w.begin(),
		 [&sum] (double x) { sum += x; return sum;});
  cout << "sum = " << sum << ", w = [ ";
  for(auto x: w) cout << x << ' '; cout << ']' << endl;
  return(0);
}
/* SAM_LISTING_END_1 */

