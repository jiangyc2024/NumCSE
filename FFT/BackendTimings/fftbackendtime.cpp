# include <complex>
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
#ifndef EIGEN_FFTW_DEFAULT
#error "fftbackendtime.cpp must be compiled with FFTW as the default backend for eigen"
#endif
// include eigen fft module kiss backend
# include "unsupported/Eigen/src/FFT/ei_kissfft_impl.h"
# include <figure/figure.hpp>
# include "timer.h"

using namespace Eigen;

//! Benchmarking different \eigen's different backends
//  for the discrete fourier transformation
//  N = maximal length of the transformed vector
//  nruns = no. of runs on which to take the minimal timing
void benchmark(const int N, const int nruns=200) {
    std::complex<double> i(0,1); // imaginary unit

    //* Version which computes for every vector size (extremely slow)
    //  MatrixXd res(N-1, 4); // save timing results and vector size
    //  for (int n = 2; n <= N; ++n) {
    //* Version which computes for powers of two only
    MatrixXd res(int(std::log2(N)), 4); // save timing results and vector size

    for (int n = 2, j = 0; n <= N; n*=2, ++j) {
        std::cout << "n: " << n << std::endl;

        // initialize input & output vector
        VectorXd y = VectorXd::Random(n);
        VectorXcd c = VectorXcd::Zero(n);

        // intialize backends
        Eigen::FFT<double, Eigen::internal::kissfft_impl<double>> fft_kiss;
        Eigen::FFT<double, Eigen::internal::fftw_impl<double>> fft_fftw;

        // warmup
        fft_kiss.fwd(c,y); 

        Timer t1; t1.start();
        for (int k = 0; k < nruns; ++k) {
            fft_kiss.fwd(c,y); 
            t1.lap();
        }

        // warmup
        fft_fftw.fwd(c,y); 

        Timer t2; t2.start();
        for (int k = 0; k < nruns; ++k) {
            fft_fftw.fwd(c,y);
            t2.lap();
        }

        // save timings
        res(j,0) = n; res(j,1) = t1.min(); 
        res(j,2) = t2.min();
        res(j,3) = t1.min()/t2.min();
    }

    mgl::Figure fig;
    fig.setlog(true, true);
    fig.plot(res.col(0), res.col(1), "m^").label("Eigen's FFT module (backend Kiss FFT)");
    fig.plot(res.col(0), res.col(2), "r<").label("Eigen's FFT module (backend FFTW)");
    fig.fplot("x*lg(x)*1e-r9","k;").label("O(n log(n))");
    fig.legend(0,1);
    fig.save("fftbackendtime");
    std::cout << res << std::endl;
}

int main() {
    // WARNING: this code may take very long to run
    benchmark(1024*8*8*8*8);
    return 0;
}
