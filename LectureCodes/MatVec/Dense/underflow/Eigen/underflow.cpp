#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <limits> 
using namespace std;
int main(){	
	cout.precision(15);
	double min = numeric_limits<double>::min();
	double res1 = M_PI*min/123456789101112;
	double res2 = res1*123456789101112/min;
	cout << res1 << endl << res2 << endl;
}
