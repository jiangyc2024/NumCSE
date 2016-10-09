#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "timer.h"
//#if INTERNAL
//#include <chrono>
//#include <figure/figure.hpp>
//#endif // INTERNAL

using namespace Eigen;

using trip = Triplet<double>;
using trip_vec = std::vector<trip>;

/* @brief Build matrix $A$ in COO format
 * @param[in] A An $n \times n$ matrix
 * @param[out] triplets The $nnz(A)$-dimensional vector of triplets forming $A$
 */
/* SAM_LISTING_BEGIN_0 */
trip_vec COO(const MatrixXd &A)
{
    // Initialization
    int m = A.rows();
	int n = A.cols();
    trip_vec triplets;

	// Number of nonzero entries
  //int nnz = A.nonzeros();
  //triplets.resize(nnz);

//#if SOLUTION
  //int k = 0;
	for(int i=0; i<m; ++i) { // Row iterator
		for(int j=0; j<n; ++j) { // Col iterator
			if(A(i,j) != 0) {
				trip triplet(i, j, A(i,j));
              //triplets[k] = triplet;
			  //++k;
				triplets.push_back(triplet);
             }
		}
	}
//#else // TEMPLATE
    // TODO: build matrix $A$ in COO format
//#endif // TEMPLATE

	return triplets;
}
/* SAM_LISTING_END_0 */

/* @brief Compute matrix product $AB$ in COO format:
 * Naive implementation.
 * @param[in] A_trips The $nnz(A)$-dimensional vector of triplets forming matrix $A$
 * @param[in] B_trips The $nnz(B)$-dimensional vector of triplets forming matrix $B$
 * @param[out] C_trips The $nnz(C)$-dimensional vector of triplets forming matrix $C = AB$
 */
/* SAM_LISTING_BEGIN_1 */
trip_vec COOprod_naive(const trip_vec &A, const trip_vec &B)
{
    // Initialization
    trip_vec C;

//#if SOLUTION
	for(auto const& a: A) {
		for(auto const& b: B) {
			if(a.col() == b.row()) {
				trip triplet(a.row(), b.col(), a.value()*b.value());
				C.push_back(triplet);
			}
		}
	}
//#else // TEMPLATE
    // TODO: compute matrix product $AB$ in COO format
//#endif // TEMPLATE

	return C;
}
// Complexity: O(nnz(A)*nnz(B))
/* SAM_LISTING_END_1 */

/* @brief Compute matrix product $AB$ in COO format:
 * Efficient implementation.
 * @param[in] A_trips The $nnz(A)$-dimensional vector of triplets forming matrix $A$
 * @param[in] B_trips The $nnz(B)$-dimensional vector of triplets forming matrix $B$
 * @param[out] C_trips The $nnz(C)$-dimensional vector of triplets forming matrix $C = AB$
 */
/* SAM_LISTING_BEGIN_2 */
trip_vec COOprod_effic(trip_vec &A, trip_vec &B)
{
	// Initialization
    trip_vec C;
	
	std::sort(A.begin(), A.end(), [&](const trip& a1, const trip& a2) {return a1.col() < a2.col();});
	std::sort(B.begin(), B.end(), [&](const trip& b1, const trip& b2) {return b1.row() < b2.row();});
	// Complexity: O(n_A*log(n_A) + n_B*log(n_B))
	
  //std::vector<int> vec(A.size() + B.size());
  //auto it = std::set_intersection(A.begin(), A.end(), B.begin(), B.end(), vec.begin(),
  //		       [&](const trip& a, const trip& b) {return a.col() < b.row();});
  //vec.resize(it - vec.begin());
	int i_A = 0, i_B = 0;
	std::set<int> intersect;
	while(i_A != A.size() && i_B != B.size())
	{
		if(A[i_A].col() < B[i_B].row()) ++i_A;
		else if(B[i_B].row() < A[i_A].col()) ++i_B;
		else {
		  intersect.insert(A[i_A].col()); // intersect.insert(B[i_B].row());
		  ++i_A; ++i_B;
		}
	}
	// Complexity: O(max(nnz(A),nnz(B)))
	
	trip_vec::iterator A_it = A.begin();
	trip_vec::iterator B_it = B.begin();
	for(auto& i: intersect) {

		A_it = std::find_if(A_it, A.end(),
					[&](const trip& a){return a.col() == i;});
		B_it = std::find_if(B_it, B.end(),
					[&](const trip& b){return b.row() == i;});
		
		while(A_it != A.end()) {
			if(A_it->col() != i) {
				break;
			} else {
				while(B_it != B.end()) {
					if(B_it->row() != i) {
						break;
					} else {
						trip triplet(A_it->row(), B_it->col(), (A_it->value())*(B_it->value()));
						C.push_back(triplet);
						++B_it;
					}
				}
				++A_it;
			}
		}
	}
	// Complexity: O(nnz(A*B))

	return C;
}
// Complexity: O(n*log(n) + O(nnz(A*B))
/* SAM_LISTING_END_2 */

int main() {
    // Initialization
    unsigned int n = 6;
    MatrixXd A(n,n), B(n,n);
    A <<  4, -1,  0, -1,  0,  0,
         -1,  4, -1,  0, -1,  0,
          0, -1,  4,  0,  0, -1,
         -1,  0,  0,  4, -1,  0,
          0, -1,  0, -1,  4, -1,
		  0,  0, -1,  0, -1,  4;
    B <<  0, -1,  2, -1,  0,  3,
         -1,  4, -1,  0, -1,  0,
          2, -1,  0,  0,  0, -1,
         -1,  0,  0,  4, -1,  0,
          0, -1,  0, -1,  0, -1,
		  3,  0, -1,  0, -1,  4;
		  
	// COO format
	trip_vec A_COO = COO(A);
	trip_vec B_COO = COO(B);

    // Compute with standard matrix multiplication and both new multipliers
    std::cout << "--> Check that the multipliers are correct" << std::endl;
    trip_vec C_eigen, C_naive, C_effic;

	C_eigen = COO(A * B);
    C_naive = COOprod_naive(A_COO, B_COO);
    C_effic = COOprod_effic(A_COO, B_COO);

  //std::cout << "Error = " << (x1 - x2).norm() << std::endl;

    // Compute runtimes of different multipliers
    std::cout << "--> Runtime comparison of naive vs efficient multiplier" << std::endl;
	// Number of repetitions
    unsigned int repeats = 3;

    // Header
    std::cout << std::setw(20) << "n"
              << std::setw(20) << "time naive [s]"
              << std::setw(20) << "time effic [s]"
              << std::endl;

    // Loop over matrix size
/*    for(unsigned int k = 4; k <= 12; ++k) {
        // Timers
        Timer tm_naive, tm_effic;
        unsigned int n = pow(2,k);

        // Repeat test many times
        for(unsigned int r = 0; r < repeats; ++r) {
            a = VectorXd::Random(n);
            b = VectorXd::Random(n);

            // Compute runtime with naive solver
            tm_naive.start();
            solveA(a,b,x);
            tm_naive.stop();
            // Compute runtime with efficient solver
            tm_effic.start();
            solveA_effic(a,b,x);
            tm_effic.stop();
        }

        // Print runtimes
        std::cout << std::setw(20) << n
                  << std::scientific << std::setprecision(3)
                  << std::setw(20) << tm_naive.min()
                  << std::setw(20) << tm_effic.min()
                  << std::endl;
    }*/
}