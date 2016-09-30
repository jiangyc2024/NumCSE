#include <iostream>

#include <vector>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using Triplet = Eigen::Triplet<double>;
using Triplets = std::vector<Triplet>;

using Vector = Eigen::VectorXd;
using Matrix = Eigen::SparseMatrix<double>;

//! \brief Solve system Ax = b with optimal complexity O(n)
//! \param[in] c Entries for matrix A
//! \param[in] b r.h.s. vector
//! \param[in] i0 index
//! \param[in] j0 index
//! \return Solution x, s.t. Ax = b
Vector solveLSE(const Vector & c, const Vector & b, unsigned int i0, unsigned int j0) {
    assert(c.size() == b.size()-1 && "Size mismatch!");
    
    // TODO: problem 6b, solve system in O(n)
}

int main(int, char**) {
    unsigned int n = 15;
    
    Triplets triplets;
    
    unsigned int ntriplets = 2*n;
    
    unsigned int i0 = 6, j0 = 4;
    
    //// PROBLEM 6a
    std::cout << "*** PROBLEM 6a:" << std::endl;
    
    // TODO: problem 6a, build sparse matrix A
}
