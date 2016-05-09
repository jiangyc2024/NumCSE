double expval(double x, double tol=1e-8){
	// Initialization
	double y = 1., term = 1., k = 1.;
	// Termination
	while(abs(term) > tol*min(y,1.)){
	term *= x/k;	// next summand
	y += term; // Summation
	++k;
	}	
	return y;	
}
