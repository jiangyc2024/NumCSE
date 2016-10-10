#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "timer.h"
// #if INTERNAL
#include <chrono>
#include <figure/figure.hpp>
// #endif // INTERNAL

using namespace Eigen;

using Trip = Triplet<double>;
using TripVec = std::vector<Trip>;

/* @brief Build matrix $A$ in COO format
 * @param[in] A An $m \times n$ matrix
 * @param[out] triplets The $nnz(A)$-dimensional vector of triplets forming $A$
 */
/* SAM_LISTING_BEGIN_0 */
TripVec Mat2COO(const MatrixXd &A)
{
    // Initialization
    int m = A.rows();
	int n = A.cols();
    TripVec triplets;

	// Number of nonzero entries
  //int nnz = A.nonzeros();
  //triplets.resize(nnz);

#if SOLUTION
  //int k = 0;
	for(int i=0; i<m; ++i) { // Row iterator
		for(int j=0; j<n; ++j) { // Col iterator
			if(A(i,j) != 0) {
				Trip triplet(i, j, A(i,j));
              //triplets[k] = triplet;
			  //++k;
				triplets.push_back(triplet);
             }
		}
	}
#else // TEMPLATE
    // TODO: build matrix $A$ in COO format
#endif // TEMPLATE

	return triplets;
}
/* SAM_LISTING_END_0 */

/* @brief Build matrix $A$ as Eigen::MatrixXd from COO format
 * @param[in] A An $nnz(A)$-dimensional vector of triplets forming matrix $A$
 * @param[out] A_mat The $m \times n$ matrix from vector of triplets $A$
 */
MatrixXd COO2Mat(const TripVec &A)
{
    // Initialization
    int m = 0, n = 0;
	for(auto const& a: A) {
		if(a.row() > m) m = a.row();
		if(a.col() > n) n = a.col();
	}
	++m; ++n; // First index is 0
	MatrixXd A_mat = MatrixXd::Zero(m, n);
	
  //A_mat.setFromTriplets(A.begin(), A.end());
	for(auto const& a: A) {
		A_mat(a.row(), a.col()) += a.value();
	}

	return A_mat;
}

/* @brief Build random binary matrix $A$ as Eigen::MatrixXd
 * @param[in] m Number of desired rows for $A$
 * @param[in] n Number of desired cols for $A$
 * @param[in] d Desired density ($nnz/(m*n)$) for $A$
 * @param[out] A An $m \times n$ random binary matrix
 */
MatrixXd randMat(int m, int n, double d)
{
    // Initialization
	MatrixXd A = MatrixXd::Zero(m, n);

	int nnz = round(m * n * d);

	for(int k=0; k<nnz; ++k) {
		int i = ceil(std::rand()/RAND_MAX*m);
		int j = ceil(std::rand()/RAND_MAX*n);

		A(i, j) = 1;
	}

	return A;
}

/* @brief Compute matrix product $AB$ in COO format:
 * Naive implementation.
 * @param[in] A The $nnz(A)$-dimensional vector of triplets forming matrix $A$
 * @param[in] B The $nnz(B)$-dimensional vector of triplets forming matrix $B$
 * @param[out] C The $nnz(C)$-dimensional vector of triplets forming matrix $C = AB$
 */
TripVec COOprod_naive(const TripVec &A, const TripVec &B)
{
    // Initialization
    TripVec C;

#if SOLUTION
	for(auto const& a: A) {
		for(auto const& b: B) {
			if(a.col() == b.row()) {
				Trip triplet(a.row(), b.col(), a.value()*b.value());
				C.push_back(triplet);
			}
		}
	}
#else // TEMPLATE
    // TODO: compute matrix product $AB$ in COO format (naive multiplier)
#endif // TEMPLATE

	return C;
}
// Complexity: O(nnz(A)*nnz(B))
/* SAM_LISTING_END_1 */

/* @brief Compute matrix product $AB$ in COO format:
 * Efficient implementation.
 * @param[in] A The $nnz(A)$-dimensional vector of triplets forming matrix $A$
 * @param[in] B The $nnz(B)$-dimensional vector of triplets forming matrix $B$
 * @param[out] C The $nnz(C)$-dimensional vector of triplets forming matrix $C = AB$
 */
/* SAM_LISTING_BEGIN_2 */
TripVec COOprod_effic(TripVec &A, TripVec &B)
{
	// Initialization
    TripVec C;
	
#if SOLUTION
	std::sort(A.begin(), A.end(), [&](const Trip& a1, const Trip& a2) {return a1.col() < a2.col();});
	std::sort(B.begin(), B.end(), [&](const Trip& b1, const Trip& b2) {return b1.row() < b2.row();});
	// Complexity: O(n_A*log(n_A) + n_B*log(n_B))
	
  //std::vector<int> vec(A.size() + B.size());
  //auto it = std::set_intersection(A.begin(), A.end(), B.begin(), B.end(), vec.begin(),
  //		      [&](const Trip& a, const Trip& b) {return a.col() < b.row();});
  //vec.resize(it - vec.begin());
	size_t i_A = 0, i_B = 0;
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
	
	TripVec::iterator A_idx = A.begin();
	TripVec::iterator B_idx = B.begin();
	for(auto i: intersect) {

		A_idx = std::find_if(A_idx, A.end(),
					[&](const Trip& a){return a.col() == i;});
		B_idx = std::find_if(B_idx, B.end(),
					[&](const Trip& b){return b.row() == i;});
		
		TripVec::iterator A_it;
		TripVec::iterator B_it;
		for(A_it=A_idx; A_it!=A.end(); ++A_it) {
			if(A_it->col() != i) {
				break;
			} else {
				for(B_it=B_idx; B_it!=B.end(); ++B_it) {
					if(B_it->row() != i) {
						break;
					} else {
						Trip triplet(A_it->row(), B_it->col(), (A_it->value())*(B_it->value()));
						C.push_back(triplet);
					}
				}
			}
		}
		A_idx = A_it;
		B_idx = B_it;
	}
#else // TEMPLATE
    // TODO: compute matrix product $AB$ in COO format (efficient multiplier)
#endif // TEMPLATE
	// Complexity: O(nnz(A*B))

	return C;
}
// Complexity: O(n*log(n) + O(nnz(A*B))
/* SAM_LISTING_END_2 */

int main() {
    // Initialization
    unsigned int n = 6;
    MatrixXd A(n,n), B(n,n);
    A << 1, 0, 0, 0, 0, 0,
         1, 0, 0, 0, 0, 0,
         1, 0, 0, 0, 0, 0,
         1, 0, 0, 0, 0, 0,
         1, 0, 0, 0, 0, 0,
		 1, 0, 0, 0, 0, 0;
    B << 1, 1, 1, 1, 1, 1,
         0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0,
		 0, 0, 0, 0, 0, 0;
		  
	// COO format
	TripVec A_COO = Mat2COO(A);
	TripVec B_COO = Mat2COO(B);

    // Compute with standard matrix multiplication and both new multipliers
    std::cout << "--> Check that the multipliers are correct" << std::endl;
    TripVec C_eigen, C_naive, C_effic;
	
	C_eigen = Mat2COO(A * B);
    C_naive = COOprod_naive(A_COO, B_COO);
    C_effic = COOprod_effic(A_COO, B_COO);

    MatrixXd C_mat_eigen;
	MatrixXd C_mat_naive;
	MatrixXd C_mat_effic;
	C_mat_eigen = COO2Mat(C_eigen);
	C_mat_naive = COO2Mat(C_naive);
	C_mat_effic = COO2Mat(C_effic);
	
	std::cout << "Error eigen vs naive = " << (C_mat_eigen - C_mat_naive).norm() << std::endl;
    std::cout << "Error naive vs effic = " << (C_mat_naive - C_mat_effic).norm() << std::endl;

    // Compute runtimes of different multipliers
    std::cout << "--> Runtime comparison of naive vs efficient multiplier" << std::endl;
	// Number of repetitions
    unsigned int repeats = 3;

// #if INTERNAL
    // sizes   will contain the size of the matrix
    // timings will contain the runtimes in seconds
    std::vector<double> sizes, timings_eigen, timings_naive, timings_effic;
// #endif
    
    // Header
    std::cout << std::setw(20) << "n"
		      << std::setw(20) << "time eigen [s]"
              << std::setw(20) << "time naive [s]"
              << std::setw(20) << "time effic [s]"
              << std::endl;

    // Seed
    int seed = static_cast<int> (std::chrono::system_clock::now().time_since_epoch().count());
    std::srand(seed);

    // Loop over matrix size
    for(unsigned int k = 4; k <= 12; ++k) {
        // Timers
        Timer tm_eigen, tm_naive, tm_effic;
        unsigned int n = pow(2,k);

        // Repeat test many times
        for(unsigned int r = 0; r < repeats; ++r) {
            // Initialization
//		    A = MatrixXd::Random(n,n);
//          B = MatrixXd::Random(n,n);
        	A = randMat(n, n, 0.1);
        	B = randMat(n, n, 0.1);
			
			// COO format
			TripVec A_COO = Mat2COO(A);
			TripVec B_COO = Mat2COO(B);

            // Compute runtime with Eigen solver
            tm_eigen.start();
            C_eigen = Mat2COO(A * B);
            tm_eigen.stop();
            // Compute runtime with naive solver
            tm_naive.start();
            C_naive = COOprod_naive(A_COO, B_COO);
            tm_naive.stop();
            // Compute runtime with efficient solver
            tm_effic.start();
            C_effic = COOprod_effic(A_COO, B_COO);
            tm_effic.stop();
        }
        
// #if INTERNAL
        // Save results in a vector
        sizes.push_back(n); // Save vector sizes
        timings_eigen.push_back(tm_eigen.min());
        timings_naive.push_back(tm_naive.min());
        timings_effic.push_back(tm_effic.min());
// #endif

        // Print runtimes
        std::cout << std::setw(20) << n
                  << std::scientific << std::setprecision(3)
				  << std::setw(20) << tm_eigen.min()
                  << std::setw(20) << tm_naive.min()
                  << std::setw(20) << tm_effic.min()
                  << std::endl;
    }
    
// #if INTERNAL
    mgl::Figure fig1;
    fig1.title("Timings of Eigen solver");
    fig1.ranges(2, 9000, 1e-8, 1e3);
    fig1.setlog(true, true); // set loglog scale
    fig1.plot(sizes, timings_eigen, " r+").label("Eigen");
    fig1.fplot("1e-9*x^3", "k|").label("O(n^3)");
    fig1.xlabel("Matrix size (n)");
    fig1.ylabel("Time [s]");
    fig1.legend(0, 1);
    fig1.save("structuredLSE_timing_eigen.eps");
    fig1.save("structuredLSE_timing_eigen.png");

    mgl::Figure fig2;
    fig2.title("Timings of naive solver");
    fig2.ranges(2, 9000, 1e-8, 1e3);
    fig2.setlog(true, true); // set loglog scale
    fig2.plot(sizes, timings_naive, " r+").label("naive");
    fig2.fplot("1e-9*x^3", "k|").label("O(n^3)");
    fig2.xlabel("Matrix size (n)");
    fig2.ylabel("Time [s]");
    fig2.legend(0, 1);
    fig2.save("structuredLSE_timing_naive.eps");
    fig2.save("structuredLSE_timing_naive.png");
    
    mgl::Figure fig3;
    fig3.title("Timings of efficient solver");
    fig3.ranges(2, 9000, 1e-8, 1e3);
    fig3.setlog(true, true); // set loglog scale
    fig3.plot(sizes, timings_effic, " r+").label("efficient");
    fig3.fplot("1e-9*x^3", "k|").label("O(n^3)");
    fig3.xlabel("Matrix size (n)");
    fig3.ylabel("Time [s]");
    fig3.legend(0, 1);
    fig3.save("structuredLSE_timing_effic.eps");
    fig3.save("structuredLSE_timing_effic.png");

    mgl::Figure fig4;
    fig4.title("Comparison of timings");
    fig4.ranges(2, 9000, 1e-8, 1e3);
    fig4.setlog(true, true); // set loglog scale
    fig4.plot(sizes, timings_eigen, " g+").label("Eigen");
    fig4.plot(sizes, timings_naive, " r+").label("naive");
    fig4.plot(sizes, timings_effic, " b+").label("efficient");
    fig4.fplot("1e-9*x^3", "k|").label("O(n^3)");
    fig4.fplot("1e-8*x^2", "k-").label("O(n^2)");
    fig4.fplot("1e-7*x", "k:").label("O(n)");
    fig4.xlabel("Matrix size (n)");
    fig4.ylabel("Time [s]");
    fig4.legend(0, 1);
    fig4.save("structuredLSE_comparison.eps");
    fig4.save("structuredLSE_comparison.png");
// #endif
}
