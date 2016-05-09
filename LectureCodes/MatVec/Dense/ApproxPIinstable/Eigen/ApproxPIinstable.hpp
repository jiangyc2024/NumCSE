//! Approximation of Pi by approximating the circumference of a
//! regular polygon
void ApproxPIinstable(double tol = 1e-8, int maxIt = 50){
	double s = sqrt(3)/2.; double An=3.*s;	// initialization (hexagon case)
	unsigned int n = 6, it = 0;				
	MatrixXd res(maxIt,4);					// matrix for storing results
	res(it,0) = n; res(it,1) = An;
	res(it,2) = An - M_PI; res(it,3)=s;
	while( it < maxIt && s > tol ){			// terminate when s is ’small enough’
		s = sqrt((1.- sqrt(1.-s*s))/2.);	// recursion for area
		n *= 2; An = n/2.*s;				// new estimate for circumference
		++it;
		res(it,0) = n; res(it,1) = An;		// store results and (absolute) error
		res(it,2) = An - M_PI; res(it,3)=s; 	
	}
	return res.topRows(it);
}
