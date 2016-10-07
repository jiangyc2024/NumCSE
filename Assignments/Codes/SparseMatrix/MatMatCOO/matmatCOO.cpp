#include <iostream>
#include <iomanip>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

/* @brief Build an $n \times n$ random matrix $A$ with density $d$ in COO format
 * \param[in] mt Random number generator
 * \param[in] n Size of square matrix $A$
 * \param[in] d Density of matrix $A$, i.e. $nnz(A)/n^2$
 * \param[out] A An $n \times n$ random matrix with density $d$ in COO format
 */
/* SAM_LISTING_BEGIN_0 */
MatrixXd buildA(std::mt19937 & mt, int n, double d);
{
	std::uniform_int_distribution<int> pos(0, n);
	std::uniform_real_distribution<double> val(-10, +10);

	int nnz = rand(d*n*n);

	// Number of nonzero entries
	for(int i=0; i<nnz; ++i) {
		pos(mt)
	}


int main() {

	// Random number generator
	std::random_device rd; // Seed
    std::mt19937 mt(rd()); // Mersenne Twister 19937 generator



}
