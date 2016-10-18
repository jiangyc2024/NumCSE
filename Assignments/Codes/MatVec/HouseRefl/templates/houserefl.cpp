#include <Eigen/Dense>
#include <iostream>
#include <ctime>

/* SAM_LISTING_BEGIN_1 */
void houserefl(const Eigen::VectorXd & v, Eigen::MatrixXd & Z)
{
    // TODO: rewrite algorithm
}
/* SAM_LISTING_END_1 */

int main(int argc, char ** argv) {
    // Check what houserefl does to random vector

    // Initialize random number generator
    srand((unsigned int) time(0));

    // Size of test vector
    unsigned int n = 6;
    // Optionally read from input arguments
    if(argc >= 2) n = std::atoi(argv[1]);

    Eigen::VectorXd v = Eigen::VectorXd::Random(n); // Not truly random if missing srand
    Eigen::MatrixXd Z;

    houserefl(v, Z);

    std::cout << "v = " << v << std::endl;
    std::cout << "Z = " << Z << std::endl;
}
