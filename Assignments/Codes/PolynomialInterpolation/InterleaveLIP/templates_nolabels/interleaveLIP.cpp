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
class PwLinIP { 
public:
	PwLinIP(const VectorXd &x, const VectorXd &t, const VectorXd &y);
	double operator()(double arg) const;
private:
	MatrixXd f;
	VectorXd tentBasCoeff(const VectorXd &x, const VectorXd &t,
						  const VectorXd &y) const;
};

/* @brief 
 * @param[in] x 
 * @param[in] t 
 * @param[in] y 
 * @param[out] s 
 */
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
	
    // TODO: 

	return s;
}

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

int main() {

}
