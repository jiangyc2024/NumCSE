//! ifference quotient approximation of the derivative of $\exp$
void diffq(){
	double h = 0.1, x = 0.0;
	for(int i = 1; i <= 16; ++i){
		double df = (std::exp(x+h) - std::exp(x))/h;
		std::cout << std::setprecision(14) << std::fixed;
		std::cout  << std::setw(5)<<  -i << std::setw(20) << df-1 << std::endl;
		h /= 10;
	}
}
