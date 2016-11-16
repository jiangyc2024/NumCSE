#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/QR>

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

void polyfit(const VectorXd &x, const VectorXd &y, VectorXd &coeff,
			 size_t order)
{
	Eigen::MatrixXd A(x.size(), order+1);
	Eigen::VectorXd result;

	assert(x.size() == y.size());
	assert(x.size() >= order + 1);

	// Create matrix
	for (size_t i=0; i<x.size(); ++i)
		for (size_t j=0; j<order+1; ++j)
			A(i, j) = pow(x(i), j);

	// Solve for linear least squares fit
	coeff = A.householderQr().solve(y);
	coeff.conservativeResize(order + 1);
}

/* @brief Compute values of interpolant in knots $\Vx$ from $(t_i,y_i)$
 * @param[in] x Vector of knots
 * @param[in] t Vector of nodes
 * @param[in] y Vector of values of interpolant in nodes $\Vt$
 * @param[out] s Vector of values of interpolant in knots $\Vx$
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd PwLineIntp(const VectorXd &x, const VectorXd &t,
					const VectorXd &y)

{
	assert(t.size() == y.size() &&
		  "t and y must have same size!");
	
	// Initialization
	size_t n = t.size();
	auto t_indices = order(t);
	size_t m = x.size();
	auto x_indices = order(x);
	// You can also implement a solution which does not need
	// sorted vectors and e.g. for each knot $x_j$ looks
	// for the closest node $t_{i1}$ and the next closest node $t_{i2}$.
	// However, such solution will not become more efficient
	// if you give as input already sorted vectors: for each knot $x_j$
	// you will always have to iterate along the sorted vector $t$
	// to find the included node $t_i$.
	
	VectorXd s(m);
	
#if SOLUTION
	size_t i = 0;
	for(size_t j=0; j<m; ++j) {
		
		bool intpOK = false;
		while(i < n-1) {
			
			bool inInterval = (t(t_indices[i]) <= x(x_indices[j])) &&
							  (x(x_indices[j]) <= t(t_indices[i+1]));
			
			if(inInterval) {
				intpOK = true;
				
				double gamma = (y(t_indices[i+1]) - y(t_indices[i])) /
							   (t(t_indices[i+1]) - t(t_indices[i]));
				double beta = y(t_indices[i]) - gamma * t(t_indices[i]);
				
				s(x_indices[j]) = gamma * x(x_indices[j]) + beta;

				break;
			} else {
				++i;
			}
		}
		if(!intpOK) {
			std::exit(EXIT_FAILURE); // $x \not\in [t_min,t_max]$
		}
	}
#else // TEMPLATE
    // TODO: 
#endif // TEMPLATE

	return s;
}
/* SAM_LISTING_END_0 */



int main() {
/* SAM_LISTING_BEGIN_1 */
// Compute convergence rate for interpolation by piecewise linear polyn.
// Uniform mesh in [0,1], singular f(t) = t^alpha, h-convergence
	
	// Initialization
	size_t NumAlph = 15;
	size_t NumN = 50;
	VectorXd alphas = VectorXd::LinSpaced(NumAlph,0.1,2.9);
	VectorXd nn = VectorXd::LinSpaced(NumN,1,50); // Used nodes
	
	// Points for evaluation and norm
	VectorXd x = VectorXd::LinSpaced(1000,0,1);
	MatrixXd s = (x.replicate(1,NumAlph)).cwiseProduct(
					alphas.transpose().replicate(x.size(),1) );

	MatrixXd Err(NumAlph,NumN); // Error with max norm
	MatrixXd LocErr(NumAlph,NumN); // Location of maximal error
	for(size_t i=0; i<NumN; ++i) {
		size_t n = nn(i);
		VectorXd t = VectorXd::LinSpaced(n+1,0,1); // Nodes
		MatrixXd y = t.replicate(1,NumAlph).cwiseProduct(
					   alphas.transpose().replicate(t.size(),1) );
				   
		for(size_t j=0; j<NumAlph; ++j) {
			VectorXd p = PwLineIntp(x, t, y.col(j)); // Interpolation
			size_t PosErr;
			Err(j,i) = (s.col(j) - p).cwiseAbs().maxCoeff(&PosErr);
			// PosErr is the index of the point in $x$ with max error
			// LocErr is the index of the subinterval with max error
			std::vector<double> tmp(t.size());
			VectorXd::Map(&tmp.front(), t.size()) = t -
						   x(PosErr)*VectorXd::Ones(t.size());
			LocErr(j,i) = count_if(tmp.begin(), tmp.end(),
							   [] (double val) {return val < 0;}) - 1;
			// Warning if the maximal error is not where expected
			if((alphas(j)<2 && LocErr(j,i)!=1) || (alphas(j)>2 && LocErr(j,i)!=n)) {
				std::cerr << "(alpha=" << alphas(j) << ", N=" << n <<
				"), max. err. in interval " << LocErr(j,i) << std::endl;
			}
		}
	}
	
std::cout << Err << std::endl;

//~ #if INTERNAL
    //~ mgl::Figure fig;
    //~ fig.title("Piecewise linear intp. on uniform meshes: error in max-norm");
    //~ fig.ranges(2, 100, 1e-5, 1);
    //~ fig.setlog(true, true); // Set loglog scale
    //~ for(size_t i=0; i<NumAlph; ++i)
		//~ fig.plot(nn, Err.row(i)).label("alpha="+alphas(j));
	//~ fig.xlabel("n = # subintervals");
    //~ fig.legend(0, 1);
    //~ fig.save("PwLineConv_cpp.eps");
//~ #endif // INTERNAL


	
/* SAM_LISTING_END_1 */



}
