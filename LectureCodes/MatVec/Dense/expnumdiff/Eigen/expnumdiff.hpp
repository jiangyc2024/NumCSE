#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include <Eigen/Dense>
#include <figure/figure.hpp>

using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
// This module provides support for multi precision floating point numbers via the MPFR C++ library which itself is built upon MPFR/GMP.
#include <unsupported/Eigen/MPRealSupport>
//! Numerical differentiation of exponential function with extended precision arithmetic
//! Uses the unsupported \eigen MPFRC++ Support module
void numericaldifferentiation(){
	// declare numeric type
	typedef mpfr::mpreal numeric_t;
	// bit number that should be evaluated
	int n = 7, k = 13;
	VectorXi bits(n); bits << 10,30,50,70,90,110,130;
	MatrixXd experr(13, bits.size()+1);
	for(int i = 0; i < n; ++i){
		// set precision to bits(i) bits (double has 53 bits)
		numeric_t::set_default_prec( bits(i) );
		numeric_t x = "1.1";// Evaluation point in extended precision \Label[line]{numdiff:1}
		for(int j=0;j<k;++j) {
			numeric_t h = mpfr::pow("2", -1-5*j);// width of difference quotient in extended precision
			// compute (absolute) error
			experr(j,i+1) = mpfr::abs(( mpfr::exp(x+h) - mpfr::exp(x) ) / h - mpfr::exp(x)).toDouble();
			experr(j,0) = h.toDouble();
		}	
	}
	// Plotting
	// ...
	/* SAM_LISTING_END_0 */
	std::cout << experr << std::endl;
	mgl::Figure fig;	fig.setlog(true, true);
	std::vector<string> col = {"b","g","r","c","m","h","q","e","u","p"};
	for(int i = 0; i < bits.size(); ++i)
		fig.plot(experr.col(0),experr.col(i+1), "+-"+col[i]).label(to_string(bits(i)) +" bits");
	//fig.title("One-sided difference quotient approx. of derivative of e^x");
	fig.xlabel("h");	fig.ylabel("error");	fig.grid(false);	fig.legend(1,0);
	fig.save("expnumdiffmultiprecision");
}

