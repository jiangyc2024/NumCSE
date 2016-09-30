#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace Eigen;

/* @brief Compute y = A*x with A banded matrix with diagonal structure
 * \param[in] a An (n-1)-dimensional vector for the upper diagonal
 * \param[in] b An (n-2)-dimensional vector for the second lower diagonal
 * \param[in] x An n-dimensional vector for Ax = y
 * \param[out] y The n-dimensional vector y = Ax
 */
/* SAM_LISTING_BEGIN_0 */
template <class Vector>
void multAx(const Vector & a, const Vector & b, const Vector & x, Vector & y) {
    unsigned int n = x.size();
    if( a.size() < n - 1 || b.size() < n - 2 ) {
        // Size check, error if do not match
        std::cerr << "Error: size mismatch!" << std::endl;
        return;
    }
    y.resize(n); // Ensure y has size n
    
#if SOLUTION
    // Handle first two rows
    if(n > 1) y(0) = 2*x(0) + a(0)*x(1);
    else { // Special case n = 1
        y(0) = 2*x(0);
        return;
    }
    if(n > 2) y(1) = 2*x(1) + a(1)*x(2);
    else { // Special case n = 2
        y(1) = 2*x(1);
        return;
    }
    
    // Row if n > 3 and without last row
    for(unsigned int i = 2; i < n-1; ++i) {
        y(i) = 2*x(i) + b(i-2)*x(i-2) + a(i)*x(i+1);
    }
    
    // Last row special case
    if(n > 2) y(n-1) = 2*x(n-1) + b(n-3)*x(n-3);
#endif // SOLUTION
}
/* SAM_LISTING_END_0 */

/* @brief Solve y = A*x with A banded matrix with upper triangular sparse structure
 * \param[in] a An (n-1)-dimensional vector for the upper diagonal
 * \param[in] r An n-dimensional vector for Ax = r
 * \param[out] x The n-dimensional vector from Ax = r
 */
/* SAM_LISTING_BEGIN_1 */
template <class Vector>
void solvelseAupper(const Vector & a, const Vector & r, Vector & x) {
    // Set up dimensions
    unsigned int n = r.size();
    x.resize(n);
    
#if SOLUTION
    // Backward substitution
    x(n-1) = 0.5 * r(n-1);
    for(int j = n-2; j >= 0; --j) {
        x(j) = 0.5*(r(j) - a(j)*x(j+1));
    }
#endif // SOLUTION
}
/* SAM_LISTING_END_1 */

/* @brief Solve y = A*x with A banded matrix using Gaussian-elimination (no pivot)
 * \param[in] a An (n-1)-dimensional vector for the upper diagonal
 * \param[in] b An (n-2)-dimensional vector for the second lower diagonal
 * \param[in] r An n-dimensional vector for Ax = r
 * \param[out] x The n-dimensional vector from Ax = r
 */
/* SAM_LISTING_BEGIN_2 */
template <class Vector>
void solvelseA(const Vector & a, const Vector & b, const Vector & r, Vector & x) {
    // Set up dimensions
    typedef typename Vector::Scalar Scalar;
    int n = r.size();
    x = r;
    
    // Fill in matrix: we reserve 5 nonzero entries per row for Gaussian fill in
    SparseMatrix<Scalar> A(n,n);
    A.reserve(5);
    for(int i = 0; i < n; ++i) {
        A.insert(i,i) = 2;
        if(i < n-1) A.insert(i,i+1) = a(i);
        if(i >= 2)  A.insert(i,i-2) = b(i-2);
    }
    A.makeCompressed();
    
#if SOLUTION
    // 1st stage: Gaussian elimination
    for(int i = 0; i < n-1; ++i) {
        for(int k = i+1; k < std::min(i+3,n); ++k) {
            Scalar fac = A.coeffRef(k,i)/A.coeffRef(i,i);
            for(int l = i; l < std::min(i+3,n); ++l) {
                A.coeffRef(k,l) -= fac * A.coeffRef(i,l);
            }
            x(k) -= fac * x(i);
        }
    }
    
    // 2nd stage: backwards substitution
    x(n-1) /= A.coeffRef(n-1,n-1);
    for(int i = n-2; i >= 0; --i) {
        for(int k = i+1; k < std::min(i+3,n); ++k) {
            x(i) -= x(k)*A.coeffRef(i,k);
        }
        x(i) /= A.coeffRef(i,i);
    }
#endif // SOLUTION
}
/* SAM_LISTING_END_2 */

/* @brief Solve y = A*x with A banded matrix using Eigen::SparseLU
 * \param[in] a An (n-1)-dimensional vector for the upper diagonal
 * \param[in] b An (n-2)-dimensional vector for the second lower diagonal
 * \param[in] r An n-dimensional vector for Ax = r
 * \param[out] x The n-dimensional vector from Ax = r
 */
/* SAM_LISTING_BEGIN_3 */
template <class Vector>
void solvelseAEigen(const Vector & a, const Vector & b, const Vector & r, Vector & x) {
    // Set up dimensions
    typedef typename Vector::Scalar Scalar;
    unsigned int n = r.size();
    
    // Fill in matrix: we reserve 3 nonzero entries per row for Gaussian fill in
    SparseMatrix<Scalar> A(n,n);
    A.reserve(3);
    for(unsigned int i = 0; i < n; ++i) {
        A.insert(i,i) = 2;
        if(i < n-1) A.insert(i,i+1) = a(i);
        if(i >= 2)  A.insert(i,i-2) = b(i-2);
    }
    A.makeCompressed();
    
#if SOLUTION
    // Call SparseLU
    SparseLU< SparseMatrix<Scalar> >   solver;
    solver.analyzePattern(A); 
    solver.compute(A);
    x = solver.solve(r);
#endif // SOLUTION
}
/* SAM_LISTING_END_3 */

int main() {
    unsigned int n = 9;
    // Compute with all three solvers
    std::cout << "*** Check that the solvers are correct" << std::endl;
    VectorXd a = VectorXd::Random(n-1);
    VectorXd b = VectorXd::Zero(n-2); // All 0s for upper diagonal structure
    VectorXd y = VectorXd::Random(n);
    VectorXd x;
    
    std::cout << "Original: " << y << std::endl;
    
    // Compute y = A*inv(A)*y
    solvelseAupper(a,y,x);
    multAx(a,b,x,y);
    
    // Should be the same as before
    std::cout << "Upper: " << y << std::endl;
    
    // Random b for own and Eigen-based solver
    b = VectorXd::Random(n-2);

    // Compute y = A*inv(A)*y
    solvelseA(a,b,y,x);
    multAx(a,b,x,y);
    
    // Should be the same as before
    std::cout << "Own: " << y << std::endl;
    
    // Compute y = A*inv(A)*y
    solvelseAEigen(a,b,y,x);
    multAx(a,b,x,y);
    
    // Should be the same as before
    std::cout << "Eigen: " << y << std::endl;
}
