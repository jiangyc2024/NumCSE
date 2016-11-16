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
	
	VectorXd s = VectorXd::Zero(m);
	
#if SOLUTION
	size_t i = 0;
	for(size_t j=0; j<m; ++j) {
		
		bool intpOK = false;
		while(i < n-1) {
			
			bool inInterval = (t(t_indices[i]) <= x(x_indices[j]) &&
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
	VectorXd s = x.replicate<1,NumAlph>().cwiseProduct(
				   alphas.transpose().replicate<x.size(),1>() );

	MatrixXd Err(NumAlph,NumN); // Error with max norm
	MatrixXd LocErr(NumAlph,NumN); // Location of maximal error
	for(size_t i=0; i<NumN; ++i) {
		size_t n = nn(i);
		VectorXd t = VectorXd::LinSpaced(n+1,0,1); // Nodes
		VectorXd y = t.replicate<1,NumAlph>().cwiseProduct(
					   alphas.transpose().replicate<t.size(),1>() );
				   
		for(size_t j=0; j<NumAlph; ++j) {
			VectorXd P = PWlineIntp(t, y.col(j), x); // Interpolation
			size_t PosErr;
			Err(j,k) = (s.col(j) - P).cwiseAbs().maxCoeff(&PosErr);
			// PosErr is the index of the point in $x$ with max error
			// LocErr is the index of the subinterval with max error
			
			bool IsOdd (int i) { return ((i%2)==1); }
			
			vector<double> tmp(t.size());
			VectorXd::Map(&tmp[0], t.size()) = t - x(PosErr);
			size_t LocErr(j,k) = count_if(tmp.begin(), tmp.end(),
								 [] (double val) {return val < 0;}) - 1;
			// Warning if the maximal error is not where expected
			if((alphas(j)<2 && LocErr(j,k)!=1) || (alphas(j)>2 && LocErr(j,k)!=n)) {
				std::cerr << "(alpha=" << alphas(j) << ", N=" << n <<
				"), max. err. in interval " << LocErr(j,k) << std::endl;
			}
		}
	}
	
	% estimate the convergence rate
rates = zeros(1,NumAlph);
for j =1:NumAlph
    ord = polyfit(log(nn), log(Err(j,:)), 1);
    % check the use of polyfit:
    %loglog([1,5], [1,5^ord(1)]*exp(ord(2))/2,'k--' );
    rates(j) = -ord(1);
end
	
	
	
figure;  
loglog(nn,Err,'.-','linewidth',2);
hold on; set(gca,'xtick',[3 5 10 20 50])
title('Piecewise linear intp. on uniform meshes: error in max-norm');
xlabel('n = # subintervals');
for j=1:NumAlph, leg{j} = sprintf('alpha=%4.2f',alphas(j)); end;
legend(leg, 'location','nwo');

% plot the convergence rate
axes('Position', [.1, .1, .3, .2]);
plot(alphas, rates,'-o','linewidth',2);
axis([alphas(1),alphas(end),0,2.2]);
xlabel('alpha');  ylabel('conv. rate');
hold on; plot([0 2],[0 2],'r', [2,alphas(end)],[2,2],'r');
print -depsc2 'PWlineConv.eps';

alpha_orders = [alphas(:),rates(:)]  % display conv. rates
	
	
	
	
/* SAM_LISTING_END_1 */



	
	
	
	
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
	
#if INTERNAL
	mgl::Figure fig;
	fig.xlabel("t");
	fig.ylabel("y");
	
	VectorXd t_left(2);
	t_left << t(0), t(1);
	VectorXd y_left(2);
	y_left << y(0), 0;
	fig.plot(t_left, y_left, "b");
	fig.plot(t_left, y_left, "b*");
	for(size_t i=1; i<n-1; ++i) {
		VectorXd t_(3);
		t_ << t(i-1), t(i), t(i+1);
		VectorXd y_(3);
		y_ << 0, y(i), 0;
		fig.plot(t_, y_, "b");
		fig.plot(t_, y_, "b*");
	}
	VectorXd t_right(2);
	t_right << t(n-2), t(n-1);
	VectorXd y_right(2);
	y_right << 0, y(n-1);
	fig.plot(t_right, y_right, "b");
	fig.plot(t_right, y_right, "b*");

	fig.title("Tent basis functions");
	fig.save("tent_basis_functions.eps");
#endif // INTERNAL
}
