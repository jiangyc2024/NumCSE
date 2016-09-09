#include <iostream>
#include <limits> 
using namespace std;
int main(){	
/* SAM_LISTING_BEGIN_0 */
cout.precision(25);
double eps = numeric_limits<double>::epsilon();
cout << fixed << 1.0 + 0.5*eps << endl
	 << 1.0 - 0.5*eps << endl
  	 << (1.0 + 2/eps) - 2/eps << endl;
/* SAM_LISTING_END_0 */
}
