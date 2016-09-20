///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>

using namespace std;

int main(){
/* SAM_LISTING_BEGIN_0 */
double x = exp(1000), y = 3/x, z = x*sin(M_PI), w = x*log(1);
cout << x << endl << y << endl << z << endl << w << endl;
/* SAM_LISTING_END_0 */
}
