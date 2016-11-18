#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace Eigen;

/*!
 * \brief order "argsort": find indices sorting a vector
 * \param values Array for which we want to find the sorting indices
 * \return Permutation of indices resulting in sorting of values
 */
std::vector<size_t> order(const VectorXd &values) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));
    std::sort(begin(indices), end(indices),
		[&](size_t a, size_t b) { return values[a] < values[b]; }
    );
    return indices;
}

/*!
 * @brief Intepolator class
 * Construct and evaluate a piecevise linear interpolation.
 */
/* SAM_LISTING_BEGIN_1 */
class PwLinIP { 
public:
    /*!
     * \brief
     * \param[in] x Vector of knots
     * \param[in] t Vector of nodes
     * \param[in] y Vector of values of interpolant in nodes
     */
	PwLinIP(const VectorXd &x, const VectorXd &t, const VectorXd &y);

    /*!
     * \brief operator() evaluate interpolant at $arg$.
     */
	double operator()(double arg) const;
private:
#if SOLUTION
	VectorXd x_;
	VectorXd t_;
	VectorXd y_;
	VectorXd s_;
#else // TEMPLATE
    // TODO: private members of intepolator class
#endif // TEMPLATE
    /*!
     * \brief Compute values of interpolant in knots $\mathbf{x}$ from $(t_i,y_i)$
     * \param[in] x Vector of knots
     * \param[in] t Vector of nodes
     * \param[in] y Vector of values of interpolant in nodes $\Vt$
     * \param[out] s Vector of values of interpolant in knots $\Vx$
     */
	VectorXd tentBasCoeff(const VectorXd &x, const VectorXd &t,
						  const VectorXd &y) const;
};
/* SAM_LISTING_END_1 */

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
	
#if SOLUTION
	// Check condition of subproblem 5.5.c
	size_t i = 0;
	size_t k = 0;
	for(size_t j=0; j<(n-1); ++j) {
		
        bool intervalOK = false;
		while(i < n) {
			
            bool inInterval = (x(x_indices[j]) < t(t_indices[i])) &&
                              (t(t_indices[i]) < x(x_indices[j+1]));
			
			if(inInterval) {
                intervalOK = true;
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
        if(!intervalOK) {
            std::cout << "I failed" << std::endl;
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
#else // TEMPLATE
    // TODO: compute interpolant in knots $\Vx$ from $(t_i,y_i)$
#endif // TEMPLATE

	return s;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_2 */
PwLinIP::PwLinIP(const VectorXd &x, const VectorXd &t,
				 const VectorXd &y)
{
	assert(t.size() == y.size() && t.size() == x.size() &&
		  "x, t, y must have same size!");
	
#if SOLUTION
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
#else // TEMPLATE
    // TODO: implement constructor of intepolator class
#endif // TEMPLATE
}
/* SAM_LISTING_END_2 */

 /* SAM_LISTING_BEGIN_4 */
double PwLinIP::operator()(double arg) const
{
#if SOLUTION
	if(arg < x_(0) || arg > x_(x_.size()-1)) {
		
		return 0;
	} else {
        std::vector<double> x(x_.size()-1); // Exclude $x_0$ from search of interval including $arg$
        VectorXd::Map(&x.front(), x_.size()-1) = x_.tail(x_.size()-1);
        size_t j = std::lower_bound(x.begin(), x.end(), arg) - x.begin() + 1; // Binary search
        // '+1' is needed to restore indexing from $x_0$

        double gamma = (s_(j) - s_(j-1)) / (x_(j) - x_(j-1));
        double beta = s_(j-1) - gamma * x_(j-1);
		
        return gamma * arg + beta;
	}
#else // TEMPLATE
	return 0; // TODO: implement operator() of intepolator class
#endif // TEMPLATE
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_3 */
int main() {
	// Initialization
	size_t n = 11;

    // Nodes
	VectorXd x = VectorXd::LinSpaced(n,0,10);
	VectorXd t(n);
	t(0) = 0; t.tail(n-1).setLinSpaced(n-1,0.5,9.5);

    mgl::Figure fig;
    fig.xlabel("t");
    fig.ylabel("y");

    // Plot colors
    const unsigned int M = 3;
    std::vector<std::string> C = {"b", "r", "m"};

    // Loop over basis functions
    for(unsigned int i = 0; i < M; ++i) {
        // Basis function $e_i$
        VectorXd y = VectorXd::Zero(n);
        y(i) = 1;

        // Create interpoland with values $y$
        PwLinIP cardinalBasis(x, t, y);

        // Plotting values
        unsigned int N = 1000;
        VectorXd xval = VectorXd::LinSpaced(N,0,10);
        VectorXd s(N);
        for(size_t j=0; j<N; ++j) {
            s(j) = cardinalBasis(xval(j));
        }

        std::stringstream ss;
        ss << "Basis k = " << i;
        fig.plot(xval, s, C[i].c_str()).label(ss.str().c_str());
        fig.legend();

        fig.title("Tent basis functions");
        fig.save("tent_basis_functions.eps");
    }

}
/* SAM_LISTING_END_3 */
