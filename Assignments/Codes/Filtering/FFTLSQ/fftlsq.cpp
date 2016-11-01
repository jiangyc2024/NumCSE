#include <iostream>

#include <Eigen/Dense>

#include "FFT/fft.hpp"

using namespace Eigen;

/*!
 * \brief gauss_fit
 * \param d
 * \param d
 * \return
 */
VectorXd gauss_fit(const VectorXd & d,
                   unsigned int m) {
    unsigned int n = d.size();

    FFT<double> fft;
    VectorXcd temp = fft.fwd(d);
    VectorXd rhs = temp.head(m).real() / n;
    rhs.tail(m-1) *= 2;

    return rhs;
}

VectorXd generate_d(unsigned int n) {

    unsigned int m = 10;

    VectorXd c = VectorXd::Random(m);

    srand(time(nullptr));
    double eps = 10e-2;
    unsigned int b = 100;

    VectorXd ret(n);
    for (unsigned int i = 0; i < n; ++i) {
        double r = 0;
        for (unsigned int j = 0; j < m; ++j) {
            r += c(j) * std::cos(2 * M_PI * i * j / n) + rand() % b / (double) b * eps;
        }
        ret(i) += r;
    }
    return ret;
}

int main(int argc, char **argv) {
    if(argc > 1) {
        std::cout << generate_d(std::stoi(argv[1]));
        return 0;
    }

    unsigned int m = 10;

    VectorXd d;
    gauss_fit(d, m);
}
