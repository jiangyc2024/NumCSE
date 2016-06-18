#include <limits>
#include <iostream>
#include <iomanip>
using namespace std;
int main(){
	cout << std::defaultfloat << numeric_limits<double>::min() << endl
	<< std::hexfloat << numeric_limits<double>::min() << endl 
	<< std::defaultfloat << numeric_limits<double>::max() << endl 
	<< std::hexfloat << numeric_limits<double>::max() << endl;
}

