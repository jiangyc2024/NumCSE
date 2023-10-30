///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <limits> 
using std::cout;
using std::endl;
int main(){	
/* SAM_LISTING_BEGIN_0 */
cout.precision(25);
const double eps = std::numeric_limits<double>::epsilon();
cout << std::fixed << 1.0 + 0.5*eps << endl
	 << 1.0 - 0.5*eps << endl
  	 << (1.0 + 2/eps) - 2/eps << endl;
/* SAM_LISTING_END_0 */
}
