///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include "compatibility.hpp" // compatibility fixes for std::defaultfloat, std::hexfloat
/* SAM_LISTING_BEGIN_0 */

#include <iomanip>
#include <iostream>
#include <limits>

using std::cout;
using std::endl;
using std::numeric_limits;

int main(){
	cout << numeric_limits<double>::is_iec559 << endl 
	<< std::defaultfloat << numeric_limits<double>::min() << endl
	<< std::hexfloat << numeric_limits<double>::min() << endl 
	<< std::defaultfloat << numeric_limits<double>::max() << endl 
	<< std::hexfloat << numeric_limits<double>::max() << endl;
}
/* SAM_LISTING_END_0 */

