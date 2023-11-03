///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

int main(){
/* SAM_LISTING_BEGIN_0 */
const double x = exp(1000);
const double y = 3/x;
const double z = x*sin(M_PI);
const double w = x*log(1);
cout << x << endl << y << endl << z << endl << w << endl;
/* SAM_LISTING_END_0 */
}
