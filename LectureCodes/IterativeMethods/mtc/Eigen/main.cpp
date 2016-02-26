#include <iostream>

double sqrtit(double a)
{
	double x_old = -1;
	double x = a;

	while (x_old != x)
	{
		x_old = x;
		x = 0.5*(x+a/x);
	}
	return x;
}

int main()
{
	std::cout << "sqrtit(2)=" << sqrtit(2) << std::endl;	
}
