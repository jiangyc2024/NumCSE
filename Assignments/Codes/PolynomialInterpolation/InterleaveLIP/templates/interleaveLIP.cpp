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
    // TODO: private members of intepolator class
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
    // sorted vectors and e.g.\ for each knot $x_j$ looks
	// for the closest node $t_{i1}$ and the next closest node $t_{i2}$.
	// However, such solution will not become more efficient
	// if you give as input already sorted vectors: for each knot $x_j$
	// you will always have to iterate along the sorted vector $t$
	// to find the included node $t_i$.
	
	VectorXd s = VectorXd::Zero(n);
	
    // TODO: compute interpolant in knots $\Vx$ from $(t_i,y_i)$

	return s;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_2 */
PwLinIP::PwLinIP(const VectorXd &x, const VectorXd &t,
				 const VectorXd &y)
{
	assert(t.size() == y.size() && t.size() == x.size() &&
		  "x, t, y must have same size!");
	
    // TODO: implement constructor of intepolator class
}
/* SAM_LISTING_END_2 */

 /* SAM_LISTING_BEGIN_4 */
double PwLinIP::operator()(double arg) const
{
	return 0; // TODO: implement operator() of intepolator class
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
