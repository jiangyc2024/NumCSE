//! C++ function computing the zeros of a quadratic polynomial
//! $\xi\to \xi^2+\alpha\xi+\beta$ by means of the familiar discriminant
//! formula $\xi_{1,2} = \frac{1}{2}(-\alpha\pm\sqrt{\alpha^2-4\beta})$. However
//! this implementation is \emph{vulnerable to round-off}! The zeros are
//! returned in a column vector
VectorXd zerosquadpol(double alpha, double beta){
	VectorXd z(2);
	double D = std::pow(alpha,2) -4*beta; // discriminant
	if(D < 0) throw "no real zeros";
	else{
		// The famous discriminant formula
		double wD = std::sqrt(D);
		z << (-alpha-wD)/2, (-alpha+wD)/2;
	}
	return z;
}
