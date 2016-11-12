#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>

using namespace Eigen;

std::vector<size_t> ordered(const VectorXd &values) {
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
 */
/* SAM_LISTING_BEGIN_1 */
class PwLinIP { 
public:
	PwLinIP(const VectorXd &x, const VectorXd &t, const VectorXd &y);
	double operator()(double arg) const;
private:
	MatrixXd f;
	VectorXd tentBasCoeff(const VectorXd &x, const VectorXd &t,
						  const VectorXd &y) const;
};
/* SAM_LISTING_END_1 */

/* @brief 
 * @param[in] x 
 * @param[in] t 
 * @param[in] y 
 * @param[out] s 
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd PwLinIP::tentBasCoeff(const VectorXd &x, const VectorXd &t,
							   const VectorXd &y) const

{
	// Initialization
	size_t m = x.size();
	size_t n = t.size();
	auto x_indices = ordered(x);
	auto t_indices = ordered(t);
	// You can also implement a solution which does not need
	// sorted vectors and e.g. for each knot $x_j$ looks
	// for the closest node $t_{i1}$ and the next closest node $t_{i2}$.
	// However, such solution will not become more efficient
	// if you give as input already sorted vectors.
	
	size_t i;
	
	i = 0;
	for(size_t j=0; j<(m-1); ++j) {
		
		bool nodeOK = false;
		while(i < n) {
			
			bool inInterval = (x(x_indices[j]) < t(t_indices[i])) &&
							  (t(t_indices[i]) < x(x_indices[j+1]));
			
			if(inInterval) {
				nodeOK = true;
				break;
			} else {
				++i;
			}
		}
		if(!nodeOK) {
			std::exit(EXIT_FAILURE);
		}
	}
	
	VectorXd s = VectorXd::Zero(n);
	
	i = 0;
	for(size_t j=0; j<m; ++j) {
		
		bool isLarger = false;
		while(i < n) {
			
			if(x(x_indices[j]) < t(t_indices[i])) {
				isLarger = true;
				break;
			} else {
				++i;
			}
			
			// Tent function of left node (descending)
			if(i != 0 && isLarger) {
				// Does not run if $x_j$ is outside of $[t_min,t_max]$
				s(x_indices[j]) += y(t_indices[i-1]) *
					   (x(x_indices[j])   - t(t_indices[i])) /
				       (t(t_indices[i-1]) - t(t_indices[i]));
			}
			// Tent function of right node (ascending)
			if(isLarger) {
				// Does not run if there is no $t_i$ which is $> x_j$
				s(x_indices[j]) += y(t_indices[i]) *
					   (x(x_indices[j]) - t(t_indices[i-1])) /
				       (t(t_indices[i]) - t(t_indices[i-1]));
			}
		}
	}

	return s;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_2 */
PwLinIP::PwLinIP(const VectorXd &x, const VectorXd &t,
				 const VectorXd &y)
{
	assert(t.size() == y.size() && "t and y must have same size!");
	
	f.resize(t.size(), 2);
	
	auto t_indices = ordered(t);
	
	for(size_t i=0; i < t.size(); ++i) {
		f.row(i) << t[t_indices[i]], y[t_indices[i]];
	}
}

double PwLinIP::operator()(double arg) const
{
	VectorXd x(1); x << arg;
	
	VectorXd s = tentBasCoeff(x, f.col(0), f.col(1));
	
	return s(0);
}
/* SAM_LISTING_END_2 */

int main() {

}
