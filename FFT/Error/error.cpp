/*
 * Calculating the error of first applying fourier transformation on a vector
 * and then inverse fourier transformation with Eigen's built-in FFT and std::vectors
 */

# include <iostream>
# include <iomanip>
# include <vector>
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>

typedef std::vector<double> real_vec;
typedef std::vector< std::complex<double> > cplx_vec;

int main()
{
  const unsigned int NN = 5e6;
  for (unsigned int N = 2; N < NN; N *= 2){
    real_vec v(N), w(N);
    cplx_vec f(N);

    for (unsigned int i = 0; i < N; ++i){
      v[i] = std::cos(5*2*M_PI*i/N);
    }

    Eigen::FFT<double> fft;
    fft.fwd(f, v);
    fft.inv(w, f);

    double err = 0;
    for (unsigned int j = 0; j < N; ++j){
      err += std::abs(w[j] - v[j]);
    }

    std::cout << "N = " << std::setw(7) << N << "  v - inv(fwd(v)) = " << err << "\n";
  }
  return 0;
}
