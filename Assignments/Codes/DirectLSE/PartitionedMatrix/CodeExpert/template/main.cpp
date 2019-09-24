
#include "blockdecomp.hpp"

int main() {
    // System dimension (-1)
    int n = 9;
    // Set random seed for reproducibility.
    std::srand(9);
    // Random test vectors
    Eigen::MatrixXd R = Eigen::MatrixXd::Random(n,n)
        .triangularView<Eigen::Upper>();
    Eigen::VectorXd v = Eigen::VectorXd::Random(n);
    Eigen::VectorXd u = Eigen::VectorXd::Random(n);
    Eigen::VectorXd bb = Eigen::VectorXd::Random(n+1);
    Eigen::VectorXd xe, xo;
    
    // Solve LSE, and test solution.
    solvelse(R, v, u, bb, xo);
    bool works = testSolveLSE(R,v,u,bb,xe);
    
    std::cout << "LSE solution: "
              << xo.transpose()
              << std::endl;
    
    if(works) {
        std::cout << "Error compared to LU: "
                  << (xe-xo).norm() << std::endl;
    } else std::cout << "testSolveLSE() returns false.";
}
