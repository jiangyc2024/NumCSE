/*
 * Plotting runtimes for fourier transformation and inverse fourier transformation
 * using Eigen's built-in FFT and Eigen::Vectors (and MathGL for plotting)
 *
 * RESULTS: for both methods (fft.fwd and fft.inv) O(n log(n)) asymptotic complexity
 */

# include <iostream>
# include <iomanip>
# include "timer.hpp"
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
# include <mgl2/mgl.h>

typedef Eigen::VectorXd real_vec;
typedef Eigen::VectorXcd cplx_vec;

int main()
{
  const unsigned int NN = 1e7;
  std::vector<double> evals, times_fwd, times_inv;
  evals.reserve(std::log(NN)/std::log(2));
  times_fwd.reserve(std::log(NN)/std::log(2));
  times_inv.reserve(std::log(NN)/std::log(2));

  for (unsigned int N = 2; N < NN; N *= 2){
    real_vec v = 100*Eigen::VectorXd::Random(N);
    cplx_vec f(N);
    Eigen::FFT<double> fft;

    Timer t_fwd, t_inv;
    t_fwd.start();
    fft.fwd(f, v);
    t_fwd.stop();

    std::cout << f(0);

    t_inv.start();
    fft.inv(v, f);
    t_inv.stop();

    std::cout << v(0);


    //std::cout << "N = " << std::setw(7) << N << " Fwd: " << std::setw(12) << t_fwd.duration() << " Inv: " << std::setw(12) << t_inv.duration() << "\n";
    evals.push_back(N);
    times_fwd.push_back(t_fwd.duration());
    times_inv.push_back(t_inv.duration());
  }

  mglData evalsd(evals.data(), evals.size()),
          fwdd(times_fwd.data(), times_fwd.size()),
          invd(times_inv.data(), times_inv.size());

  mglGraph gr;
  gr.SetFontSizePT(8);
  gr.SetRanges(2, NN, invd.Minimal(), fwdd.Maximal());
  gr.SetFunc("lg(x)", "lg(y)");
  gr.Axis();
  gr.Grid("","h");
  gr.Label('x', "Vector length n", 0);
  gr.Label('y', "Runtime t in s", 0);

  gr.Plot(evalsd, fwdd, "r0");
  gr.Plot(evalsd, invd, "b0");
  gr.FPlot("x*lg(x)/1e6","k;");

  gr.AddLegend("fft(v)", "r0");
  gr.AddLegend("ifft(f)", "b0");
  gr.AddLegend("O(n log(n))", "k;");
  gr.Legend(1,0);

  gr.WriteEPS("runtimes.eps");

  return 0;
}

    

