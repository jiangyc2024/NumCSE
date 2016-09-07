#include <iostream>
int main(){
	std::cout.precision(15);
	double a = 4.0/3.0, b = a-1, c = 3*b, e = 1-c;
	std::cout << e << std::endl;
	a = 1012.0/113.0; b = a-9; c = 113*b; e = 5+c;
	std::cout << e << std::endl;
	a = 83810206.0/6789.0; b = a-12345; c = 6789*b; e = c-1;
	std::cout << e << std::endl;
}

