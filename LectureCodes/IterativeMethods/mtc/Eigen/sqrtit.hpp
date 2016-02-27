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
