#pragma once
#pragma begin<0>
// This module provides support for multi precision floating point numbers
// via the MPFR C++ library which itself is built upon MPFR/GMP.
#include <unsupported/Eigen/MPRealSupport>
//! Numerical differentiation of exponential function with extended precision arithmetic
//! Uses the unsupported Eigen MPFRC++ Support module
void numericaldifferentiation(){
	typedef mpfr::mpreal numeric_t;
	// Declare matrix and vector types with multi-precision scalar type
	typedef Matrix<numeric_t,Dynamic,Dynamic>  MatrixXmp;
	typedef Matrix<numeric_t,Dynamic,1>        VectorXmp;
	// bit number that should be evaluated
	int n = 7, k = 13;
	VectorXi bits(n); bits << 10,30,50,70,90,110,130;
	MatrixXd experr(13, bits.size()+1);
	for(int i = 0; i < n; ++i){
		// set precision to bits(i) bits (double has 53 bits)
		mpfr::mpreal::set_default_prec( bits(i) );
		mpfr::mpreal x = "1.1";// Evaluation point in extended precision
		for(int j=0;j<k;++j) {
			mpfr::mpreal h = mpfr::pow("2", -1-5*j);// width of difference quotient in extended precision
			// compute (absolute) error
			experr(j,i+1) = mpfr::abs(( mpfr::exp(x+h) - mpfr::exp(x) ) / h - mpfr::exp(x)).toDouble();
			experr(j,0) = h.toDouble();
		}	
	}
	// Plotting
	mgl::Figure fig;	fig.setFontSize(4);	fig.setlog(true, true);
	std::vector<string> col = {"b","g","r","c","m","h","q","e","u","p"};
	for(int i = 0; i < bits.size(); ++i)
		fig.plot(experr.col(0),experr.col(i+1), "+-"+col[i]).label(to_string(bits(i)) +" bits");
	fig.title("One-sided difference quotient approx. of derivative of e^x");
	fig.xlabel("h");	fig.ylabel("error");	fig.grid(false);	fig.legend(1,0);
	fig.save("expnumdiffmultiprecision");
}
#pragma end<0>
