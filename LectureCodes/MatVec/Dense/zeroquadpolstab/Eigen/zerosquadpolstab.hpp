//! C++ function computing the zeros of a quadratic polynomial
//! $\xi\to \xi^2+\alpha\xi+\beta$ by means of the familiar discriminant
//! formula $\xi_{1,2} = \frac{1}{2}(-\alpha\pm\sqrt{\alpha^2-4\beta})$. 
//! This is a stable implementation based on Vieta's theorem.
//! The zeros are returned in a column vector
VectorXd zerosquadpolstab(double alpha, double beta){
	VectorXd z(2);
	double D = std::pow(alpha,2) -4*beta; // discriminant
	if(D < 0) throw "no real zeros";
	else{
		double wD = std::sqrt(D);
		// Use discriminant formula only for zero far away from $0$ 
		// in order to \com{avoid cancellation}. For the other zero
		// use Vieta's formula. 
		if(alpha >= 0){
			double t = 0.5*(-alpha-wD);	// \Label[line]{zqs:11}
			z << t, beta/t;			
		}
		else{
				double t = 0.5*(-alpha+wD);	// \Label[line]{zqs:12}
				z << beta/t, t;
		}
	}
	return z;
}
