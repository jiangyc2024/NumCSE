//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace Eigen;

std::vector<size_t> order(const VectorXd &values) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));
    std::sort(begin(indices), end(indices),
		[&](size_t a, size_t b) { return values[a] < values[b]; }
    );
    return indices;
}

/* @brief Intepolator class
 * @param[in] x Vector of knots
 * @param[in] t Vector of nodes
 * @param[in] y Vector of values of interpolant in nodes
 */
class PwLinIP { 
public:
	PwLinIP(const VectorXd &x, const VectorXd &t, const VectorXd &y);
	double operator()(double arg) const;
private:
	VectorXd x_;
	VectorXd t_;
	VectorXd y_;
	VectorXd s_;
	VectorXd tentBasCoeff(const VectorXd &x, const VectorXd &t,
						  const VectorXd &y) const;
};

/* @brief Compute values of interpolant in knots $\Vx$ from $(t_i,y_i)$
 * @param[in] x Vector of knots
 * @param[in] t Vector of nodes
 * @param[in] y Vector of values of interpolant in nodes $\Vt$
 * @param[out] s Vector of values of interpolant in knots $\Vx$
 */
VectorXd PwLinIP::tentBasCoeff(const VectorXd &x, const VectorXd &t,
							   const VectorXd &y) const

{
	// Initialization
	size_t n = t.size();
	auto x_indices = order(x);
	auto t_indices = order(t);
	// You can also implement a solution which does not need
	// sorted vectors and e.g. for each knot $x_j$ looks
	// for the closest node $t_{i1}$ and the next closest node $t_{i2}$.
	// However, such solution will not become more efficient
	// if you give as input already sorted vectors: for each knot $x_j$
	// you will always have to iterate along the sorted vector $t$
	// to find the included node $t_i$.
	
	VectorXd s = VectorXd::Zero(n);
	
    // TODO: 

	return s;
}

PwLinIP::PwLinIP(const VectorXd &x, const VectorXd &t,
				 const VectorXd &y)
{
	assert(t.size() == y.size() && t.size() == x.size() &&
		  "x, t, y must have same size!");
	
	size_t n = t.size();
	x_.resize(n);
	t_.resize(n);
	y_.resize(n);
	
	auto x_indices = order(x);
	for(size_t i=0; i<n; ++i) {
		x_(i) = x[x_indices[i]];
	}
	
	auto t_indices = order(t);
	for(size_t i=0; i<n; ++i) {
		t_(i) = t[t_indices[i]];
		y_(i) = y[t_indices[i]];
	}
	
	s_ = tentBasCoeff(x_, t_, y_);
}

double PwLinIP::operator()(double arg) const
{	
	if(arg < x_(0) || arg > x_(x_.size()-1)) {
		
		return 0;
	} else {
		
		size_t j = 1; // Already checked that arg >= x_(0)
		while(j < x_.size()) {
			if(arg <= x_(j)) {
				break;
			} else {
				++j;
			}
		}
		
		double gamma = (s_(j) - s_(j-1)) / (x_(j) - x_(j-1));
		double beta = s_(j-1) - gamma * x_(j-1);

		return gamma * arg + beta;
	}
}

int main() {
	// Initialization
	size_t n = 11;
	VectorXd x = VectorXd::LinSpaced(n,0,10);
	VectorXd t(n);
	t(0) = 0; t.tail(n-1).setLinSpaced(n-1,0.5,9.5);
	VectorXd y = VectorXd::Ones(n);
	
	PwLinIP cardinalBasis(x, t, y);
	
	VectorXd s(n);
	for(size_t j=0; j<n; ++j) {
		s(j) = cardinalBasis(x(j));
	}
	
}
