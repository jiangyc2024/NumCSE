#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

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
    int n = A.rows();
    trip_vec triplets;

	// Number of nonzero entries
	int nnz = A.nonzeros();
	triplets.resize(nnz);

#if SOLUTION
	int k = 0;
	for(int i=0; i<m; ++i) { // Row iterator
		for(int j=0; j<n; ++j) { // Col iterator
			if(A(i,j) != 0) {
				trip triplet(i, j, A(i,j));
              	triplets(k) = triplet;
				++k;
             }
		}
	}
#else // TEMPLATE
    // TODO: build matrix $A$ in COO format
#endif // TEMPLATE

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
trip_vec COO_prod(const trip_vec &A, const trip_vec &B)
{
    // Initialization
    trip_vec C;

	for(auto const& a: A) {
		for(auto const& b: B) {
			if(a.col() == b.row()) {
				trip triplet(a.row(), b.col(), a.value*b.value());
				C.push_back(triplet);
			}
		}
     }

	return C;
}
// Complexity: O(nnz(A)*nnz(B))
/* SAM_LISTING_END_1 */

/* @brief Compute matrix product $AB$ in COO format:
 * Naive implementation.
 * @param[in] A_trips The $nnz(A)$-dimensional vector of triplets forming matrix $A$
 * @param[in] B_trips The $nnz(B)$-dimensional vector of triplets forming matrix $B$
 * @param[out] C_trips The $nnz(C)$-dimensional vector of triplets forming matrix $C = AB$
 */
/* SAM_LISTING_BEGIN_2 */
trip_vec COO_prod(const trip_vec &A, const trip_vec &B)
{
	std::sort(A.begin(), A.end(), [&](const trip& a1, const trip& a2) {return a1.col() < a2.col();});
	std::sort(B.begin(), B.end(), [&](const trip& b1, const trip& b2) {return b1.row() < b2.row();});
	// Complexity: O(n_A*log(n_A) + n_B*log(n_B))

	std::vector<int> vec(A.size() + B.size());
  	auto it = std::set_intersection(A.begin(), A.end(), B.begin(), B.end(), vec.begin(),
			      [&](const auto& a, const auto& b) {return a.col() < b.row();});
	vec.resize(it - vec.begin());

	for(it=vec.begin(); it!=vec.end(); ++it) {
		
		auto A_it = std::find_if(std::begin(A), std::end(A),
						[&](const auto& a){return a.col() == *it;});
		auto B_it = std::find_if(std::begin(B), std::end(B),
						[&](const auto& b){return b.row() == *it;});

		for(auto& a_it: A_it) {
			if(*a_it.col() != *it) {
				break;
			}
			for(auto& b_it: B_it) {
				if(*b_it.row() != *it) {
					break;
				} else {
					trip triplet(*a_it.row(), *b_it.col(), (*a_it.value)*(*b_it.value()));
				C.push_back(triplet);
			}
		}
	}
	// Complexity: O(nnz(A*B))

	return C;
}
// Complexity: O(n*log(n) + O(nnz(A*B))
/* SAM_LISTING_END_2 */

