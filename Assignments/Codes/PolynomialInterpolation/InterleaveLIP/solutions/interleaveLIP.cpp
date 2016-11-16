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
/* SAM_LISTING_BEGIN_1 */
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
/* SAM_LISTING_END_1 */

/* @brief Compute values of interpolant in knots $\Vx$ from $(t_i,y_i)$
 * @param[in] x Vector of knots
 * @param[in] t Vector of nodes
 * @param[in] y Vector of values of interpolant in nodes $\Vt$
 * @param[out] s Vector of values of interpolant in knots $\Vx$
 */
/* SAM_LISTING_BEGIN_0 */
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
	
	// Check condition of subproblem 5.5.c
	size_t i = 0;
	size_t k = 0;
	for(size_t j=0; j<(n-1); ++j) {
		
		bool nodeOK = false;
		while(i < n) {
			
			bool inInterval = (x(x_indices[j]) <= t(t_indices[i])) &&
							  (t(t_indices[i]) <= x(x_indices[j+1]));
			
			if(inInterval) {
				nodeOK = true;
				if(i == j) { // Index of interval which contains 2 nodes
							 // $t_i$ and $t_{i+1}$. After that, we have
							 // $i > j$...
					k = j;
				}
				break;
			} else {
				++i;
			}
		}
		if(!nodeOK) {
			std::exit(EXIT_FAILURE);
		}
	}
	
	// 1. Find slope $\gamma$ and intercept $\beta$
	// in interval with 2 nodes $k$
	// 2. Find $s_k$ and $s_{k+1}$
	double gamma = (y(t_indices[k+1]) - y(t_indices[k])) /
				   (t(t_indices[k+1]) - t(t_indices[k]));
	double beta = y(t_indices[k]) - gamma * t(t_indices[k]);
	
	s(x_indices[k])   = gamma * x(x_indices[k])   + beta;
	s(x_indices[k+1]) = gamma * x(x_indices[k+1]) + beta;
	
	// Find intercept, slope and value $s$ at the lower bound $x$
	// for all intervals on the left of interval with 2 nodes $k$
	for(int j=k-1; j>=0; --j) {
		gamma = (s(x_indices[j+1]) - y(t_indices[j])) /
				(x(x_indices[j+1]) - t(t_indices[j]));
		beta = y(t_indices[j]) - gamma * t(t_indices[j]);
		
		s(x_indices[j]) = gamma * x(x_indices[j]) + beta;
	}
	// Find intercept, slope and value $s$ at the upper bound $x$
	// for all intervals on the right of interval with 2 nodes $k$
	for(int j=k+2; j<n; ++j) {
		gamma = (y(t_indices[j]) - s(x_indices[j-1])) /
				(t(t_indices[j]) - x(x_indices[j-1]));
		beta = s(x_indices[j-1]) - gamma * x(x_indices[j-1]);
		
		s(x_indices[j]) = gamma * x(x_indices[j]) + beta;
	}

	return s;
}
/* SAM_LISTING_END_0 */

/* @brief Constructor of intepolator class
 */
/* SAM_LISTING_BEGIN_2 */
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
/* SAM_LISTING_END_2 */

/* @brief Operator() of intepolator class
 */
 /* SAM_LISTING_BEGIN_4 */
double PwLinIP::operator()(double arg) const
{
	if(arg < x_(0) || arg > x_(x_.size()-1)) {
		
		return 0;
	} else {
		size_t j = 1; // Already checked that $arg \geq x_0$
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
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_3 */
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
/* SAM_LISTING_END_3 */
