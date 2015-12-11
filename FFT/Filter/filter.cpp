# include <iostream>
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
# include <mgl2/mgl.h>

void plot(mglGraph* gr, const Eigen::VectorXd& t, const Eigen::VectorXd& y, const char* style, const char* label)
{
  mglData yd(y.data(), y.size()),
          td(t.data(), t.size());

  gr->SetRanges(td.Minimal(), td.Maximal()/2., yd.Minimal(), yd.Maximal());
  gr->SetOrigin(mglPoint(0,0));
  gr->Axis();
  gr->Plot(td, yd, style);

  gr->AddLegend(label, style);
  gr->Legend();
  return;
}

void filter(Eigen::VectorXd& signal, Eigen::VectorXd& model)
{
  const long N = signal.size();
  Eigen::FFT<double> fft;
  Eigen::VectorXcd k(N), c(N);
  // transform to spectrum of frequencies
  fft.fwd(k, signal);
  // get strong frequencies
  const double T = k.cwiseAbs().maxCoeff()/2.; // threshold
  for (long n = 0; n < N; ++n){
    if (std::abs(k(n)) > T){
      c(n) = k(n);
    }
    else {
      c(n) = 0.;
    }
  }
  // transform back with inverse fourier transform
  fft.inv(model, c);
  Eigen::VectorXd f = Eigen::VectorXd::LinSpaced(N, 1, N);

  mglGraph gr;
  gr.SetFontSizePT(7);
  gr.SetRanges(0, N/2, 0, c.cwiseAbs().maxCoeff());
  gr.Title("Spectrum of frequencies");
  gr.SubPlot(2,1,0,"<_");

  Eigen::VectorXd k_abs = k.cwiseAbs(); // cant use k.cwiseAbs().data() need
  Eigen::VectorXd c_abs = c.cwiseAbs(); // a copy of k.cwiseAbs() to call data()

  mglData fd(f.data(), f.size()),
          kd(k_abs.data(), k.size()),
          cd(c_abs.data(), c.size());

  gr.Grid("","h");
  gr.Axis();
  gr.Label('x', "Frequencies f",0);
  gr.Label('y', "Amplitudes |k(f)|",0);

  gr.Plot(fd, kd, "r0");
  gr.AddLegend("Distorted", "r0");
  gr.Legend(1,0.8);

  gr.SubPlot(2,1,1,"<_");
  gr.Grid("","h");
  gr.Axis();
  gr.Label('x', "Frequencies f",0);
  gr.Label('y', "Amplitudes |k(f)|",0);
  gr.Plot(fd, cd, "b0");
  gr.ClearLegend();
  gr.AddLegend("Clean", "b0");
  gr.Legend(1,0.8);

  gr.WriteEPS("freqs.eps");

  return;
}

int main()
{
  const long N = 1000;
  const unsigned int f1(50), f2(120); // set two frequencies
  const double a1(1), a2(0.7); // set two amplitudes
  Eigen::VectorXd signal(N), dist_signal(N), model(N), t(N), noise(N);
  t = Eigen::VectorXd::LinSpaced(N, 0, 1 - 1./N);
  noise = Eigen::VectorXd::Random(N);
  // signal = a1*sin(2pi*f1*t) + a2*sin(2pi*f2*t)
  signal = ( a1*(f1*2*M_PI*t.array()).sin() + a2*(f2*2*M_PI*t.array()).sin() ).matrix();
  // dist_signal = signal + noise
  dist_signal = (signal.array() + 5*noise.array()).matrix();
  filter(dist_signal, model);

  // prepare data for plotting
  mglData td(t.data(), t.size()),
          signald(signal.data(), signal.size()),
          dist_signald(dist_signal.data(), dist_signal.size()),
          modeld(model.data(), model.size());

  mglGraph gr;
  gr.SetFontSizePT(7);

  // plot distorted signal and signal
  gr.SubPlot(1,2,0,"<_");

  gr.SetRanges(0, 0.5, -7, 7);
  gr.Grid("","h");
  gr.Axis();
  gr.Label('x', "Time t",0);
  gr.Label('y', "Signal f(t)",0);

  gr.Plot(td, dist_signald, "r0");
  gr.Plot(td, signald, "b0");

  gr.AddLegend("Original signal", "b0");
  gr.AddLegend("Distorted signal", "r0");
  gr.Legend();

  // plot signal and model
  gr.SubPlot(1,2,1,"<_");
  gr.SetRanges(0,0.2,-2.3,2.3);
  gr.Grid("","h");
  gr.Axis();
  gr.Label('x', "Time t",0);
  gr.Label('y', "Signal f(t)",0);

  gr.Plot(td, signald, "b0");
  gr.Plot(td, modeld, "g0");

  gr.ClearLegend();
  gr.AddLegend("Original signal", "b0");
  gr.AddLegend("Model", "g0");
  gr.Legend();

  gr.WriteEPS("signal.eps");

  return 0;
}
