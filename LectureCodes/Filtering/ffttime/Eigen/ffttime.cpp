# include <complex>
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
# include <figure/figure.hpp>
# include "timer.h"
# include "meshgrid.hpp"
# include "eigen_fft_backend_name.hpp"

using namespace Eigen;

//! Benchmarking discrete fourier transformations
//  N = maximal length of the transformed vector
//  nruns = no. of runs on which to take the minimal timing
//  If the size of the transformed array is a prime number the kissfft
//  implementation of \eigen{} performs extremly bad, use FFTW or MKL
//  in these cases
void benchmark(const int N, const int nruns=3) {
  std::complex<double> i(0,1); // imaginary unit

//* Version which computes for every vector size (extremely slow)
//  MatrixXd res(N-1, 4); // save timing results and vector size
//  for (int n = 2; n <= N; ++n) {
//* Version which computes for powers of two only
  MatrixXd res(int(std::log2(N)), 4); // save timing results and vector size
  for (int n = 2, j = 0; n <= N; n*=2, ++j) {
    VectorXd y = VectorXd::Random(n);
    VectorXcd c = VectorXcd::Zero(n);
    // compute fourier transformed with a loop implementation
    Timer t1; t1.start();
    for (int k = 0; k < nruns; ++k) {
      std::complex<double> omega = std::exp(-2*M_PI/n*i), 
                               s = omega;
      c(0) = y.sum();
      for (int j = 1; j < n; ++j) {
        c(j) = y(n-1);
        for (int l = n-1; l > 0; --l){ c(j) = c(j)*s + y(l); }
        s *= omega;
      }
      t1.lap();
    }
    // compute fourier transformed through matrix multiplication
    MatrixXd I, J; VectorXd lin = VectorXd::LinSpaced(n, 0, n-1);
    meshgrid(lin, lin, I, J);
    MatrixXcd F = (-2*M_PI/n*i*I.cwiseProduct(J)).array().exp().matrix();
    Timer t2; t2.start();
    for (int k = 0; k < nruns; ++k) { 
      c = F*y; t2.lap();
    }
    // use Eigen's FFT to compute fourier transformation
    // Note: slow for large primes!
    Eigen::FFT<double> fft;
    Timer t3; t3.start();
    for (int k = 0; k < nruns; ++k) { 
      fft.fwd(c,y); t3.lap();
    }
    // save timings
    res(j,0) = n; res(j,1) = t1.min(); 
    res(j,2) = t2.min(); res(j,3) = t3.min();
  }

  mgl::Figure fig;
  fig.title("FFT timing");
  fig.setlog(true, true);
  fig.plot(res.col(0), res.col(1), "b*").label("Loop based computation");
  fig.plot(res.col(0), res.col(2), "m^").label("Direct matrix multiplication");
  fig.plot(res.col(0), res.col(3), "r<").label("Eigen's FFT module");
  fig.fplot("x^2*1e-8","k:").label("O(n^2)");
  fig.fplot("x*lg(x)*1e-8","k;").label("O(n log(n))");
  fig.legend(0,1);
  fig.save("ffttime");
  std::cout << res << std::endl;
}

int main() {
  // print the backend used by eigen for the fft
  std::cout << "FFT Backend: " << eigen_fft_backend_name() << std::endl;
    
  // WARNING: this code may take very long to run
  benchmark(1024*8);
  return 0;
}
