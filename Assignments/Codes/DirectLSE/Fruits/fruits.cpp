#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

using namespace Eigen;

/* @brief Solve the LSE arising from
 * the problem description to determine the
 * prices of the fruits using Eigen
 */
/* SAM_LISTING_BEGIN_1 */
int main() {
#if SOLUTION
    // Initialize the matrix
    MatrixXd A(6,6);
    // See "Advanced Initialization" in the Eigen doc
    A << 3, 1, 7, 2, 0, 0,
         2, 0, 0, 0, 0, 1,
         1, 0, 0, 0, 3, 0,
         0, 5, 0, 1, 0, 0,
         0, 0, 0, 2, 0, 1,
         0, 1,20, 0, 0, 0;
    
    // Initialize the vector
    VectorXd b(6);
    b << 11.10, 17.00, 6.10, 5.25, 12.50, 7.00;
    
    // Initialize the solver
    ColPivHouseholderQR<MatrixXd> solver;
    // Decompose A
    solver.compute(A);
    // Solve the LSE
    VectorXd x1 = solver.solve(b);
    
    // Notice how we split up initialization, decompositon and solving.
    // This means we can use the same solver for multiple matrices and the
    // same decompositions for multiple LSE with the same matrix.
    // These three steps can be combined:
    VectorXd x2 = A.colPivHouseholderQr().solve(b);
    
    // The result is the same:
    assert((x1-x2).norm() < 1e-10);
    
    //Print the result nicely
    std::vector<std::string> fruit_names {"Mango", "Kiwi", "Lychee", "Banana", "Pomegranate", "Pineapple"};
    std::cout << std::fixed << std::setprecision(2);    //stream formatting
    std::cout <<"The fruits cost:" << std::endl;
    for (int i = 0; i < 6; i++)
    {
        std::cout << fruit_names[i] << ": " << x1[i] << " sFr." << std::endl;
    }
#else
    // Initialize an appropriate Matrix A and Vector b
    // and solve the system Ax=b
#endif
}
/* SAM_LISTING_END_1 */
