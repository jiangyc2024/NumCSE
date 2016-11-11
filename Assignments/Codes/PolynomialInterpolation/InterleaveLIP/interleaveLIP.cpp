#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

#include "timer.h"
#if INTERNAL
#include <chrono>
#include <figure/figure.hpp>
#endif // INTERNAL

using namespace Eigen;

template <typename T>
std::vector<size_t> ordered(std::vector<T> const& values) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));
    std::sort(begin(indices), end(indices),
		[&](size_t a, size_t b) { return values[a] < values[b]; }
    );
    return indices;
}

/* @brief 
 * @param[in] x
 * @param[in] t
 * @param[in] y
 * @param[out] s
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd tentBasCoeff(const VectorXd &x, const VectorXd &t,
					  const VectorXd &y);

{
	assert(t.size() == y.size() && "t and y must have same size!")
	
	// Initialization
	size_t n = x.size();
	size_t m = t.size();
	auto x_indices = ordered(x);
	auto t_indices = ordered(t);
	
#if SOLUTION
	size_t j;
	
	j = 0;
	for(size_t i=0; i<(n-1); ++i) {
		
		bool nodeOK = false;
		while(j < m) {
			
			bool inInterval = (x(x_indices[i]) < t(t_indices[j])) &&
							  (t(t_indices[j]) < x(x_indices[i+1]));
			++j;
			
			if(inInterval) {
				nodeOK = true;
				break;
			}
		}
		if(!nodeOK) {
			std::exit(EXIT_FAILURE);
		}
	}
	
	VectorXd s = VectorXd::Zero(n);
	
	j = 0;
	for(size_t i=0; i<n; ++i) {
		
		bool isLarger = false;
		while(j < m) {
			
			if(x(x_indices[i]) < t(t_indices[j])) {
				isLarger = true;
				break;
			} else {
				++j;
			}
			
			// Tent function of left node (descending)
			if(j != 0 && isLarger) {
				// Does not work if $x_i$ is outside of $[t_min,t_max]$
				s(x_indices[i]) += y(t_indices[j-1]) *
					   (x(x_indices[i])   - t(t_indices[j])) /
				       (t(t_indices[j-1]) - t(t_indices[j]));
			}
			// Tent function of right node (ascending)
			if(isLarger) {
				// Does not work if there is no $t_j$ which is $> x_i$
				s(x_indices[i]) += y(t_indices[j]) *
					   (x(x_indices[i]) - t(t_indices[j-1])) /
				       (t(t_indices[j]) - t(t_indices[j-1]));
			}
		}
	}
#else // TEMPLATE
    // TODO: 
#endif // TEMPLATE

	return s;
}
/* SAM_LISTING_END_0 */



int main() {
	// Initialization
	int m = 4;
	int n = 3;
	VectorXd u(m);
	VectorXd v(n);
	u << 1, 2, 3, 4;
	v << 10, 20, 30;

	// Compute with both functions
	std::cout << "Check that all functions are correct" << std::endl;

	VectorXd uv_1 = polyMult_naive(u, v);
	std::cout << "Naive multiplicator: "
			  << std::endl << uv_1 << std::endl;

	VectorXd uv_2 = polyMult_fast(u, v);
	std::cout << "Efficient multiplicator: "
			  << std::endl << uv_2 << std::endl;

	std::cout << "Error = " << (uv_1 - uv_2).norm() << std::endl;

	VectorXd v_new = polyDiv(uv_2, u);
	std::cout << "Error of efficient division = " << (v - v_new).norm()
			  << std::endl;

	// Initialization
	int repeats = 3;
	VectorXd uv;
#if INTERNAL
    // sizes   will contain the size of the vectors
    // timings will contain the runtimes in seconds
    std::vector<double> sizes, timings_naive, timings_effic;
#endif

	// Compute runtimes of different multiplicators
	std::cout << "Runtime comparison of naive v efficient multiplicator"
			  << std::endl;

	// Header
	std::cout << std::setw(20) << "n"
			  << std::setw(20) << "time naive [s]"
			  << std::setw(20) << "time fast [s]"
			  << std::endl;

	// Loop over vector size
	for(int p = 2; p <= 15; ++p) {
		// Timers
		Timer tm_naive, tm_effic;
		int n = pow(2,p);

		// Repeat test many times
		for(int r = 0; r < repeats; ++r) {
			u = VectorXd::Random(n);
			v = VectorXd::Random(n);

			// Compute runtime of naive multiplicator
			tm_naive.start();
			uv = polyMult_naive(u, v);
			tm_naive.stop();
			// Compute runtime of efficient multiplicator
			tm_effic.start();
			uv = polyMult_fast(u, v);
			tm_effic.stop();
		}
		
#if INTERNAL
        // Save results in a vector
        sizes.push_back(n); // Save vector sizes
        timings_naive.push_back(tm_naive.min());
        timings_effic.push_back(tm_effic.min());
#endif
		
		// Print runtimes
		std::cout << std::setw(20) << n
				  << std::scientific << std::setprecision(3)
				  << std::setw(20) << tm_naive.min()
				  << std::setw(20) << tm_effic.min()
				  << std::endl;
	}
	
#if INTERNAL
    mgl::Figure fig;
    fig.title("Comparison of timings of polynomial multiplication");
    fig.ranges(1, 35000, 1e-8, 1e3);
    fig.setlog(true, true); // Set loglog scale
    fig.plot(sizes, timings_naive, " r+").label("naive");
    fig.plot(sizes, timings_effic, " b+").label("efficient");
    fig.fplot("1e-9*x^3", "k|").label("O(n^3)");
    fig.fplot("1e-8*x^2", "k-").label("O(n^2)");
    fig.fplot("1e-7*x", "k:").label("O(n)");
    fig.xlabel("Vector size (n)");
    fig.ylabel("Time [s]");
    fig.legend(0, 1);
    fig.save("polydiv_comparison.eps");
    fig.save("polydiv_comparison.png");
#endif
}
