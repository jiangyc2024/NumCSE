# include <Eigen/Dense>
# include <mgl2/mgl.h>
# undef I
# include <chrono>


typedef std::chrono::high_resolution_clock hrclock;
double sample(const unsigned int n, const bool write)
{
  Eigen::VectorXd init = Eigen::VectorXd::Random(n);
  mglGraph gr;
  
  auto t = hrclock::now();
  mglData initd(init.data(), init.size());
  gr.Plot(initd);
  if (write)
    gr.WriteEPS("buffer.eps");
  return std::chrono::duration_cast<std::chrono::nanoseconds>(hrclock::now() - t).count();
}

void measure()
{
  const unsigned int N = 1e6;
  std::vector<double> evals; 
  std::vector<double> t_plot;
  std::vector<double> t_plotwrite;
  t_plot.reserve(log(1e6)/log(2));
  t_plotwrite.reserve(log(1e6)/log(2));
  evals.reserve(log(1e6)/log(2));

  for (unsigned int n = 2; n < N; n *= 2){
    evals.push_back(n);
    t_plot.push_back(sample(n, false));
    t_plotwrite.push_back(sample(n, true));
  }
  
  mglGraph gr;
  mglData evalsd(evals.data(), evals.size());
  mglData t_plotd(t_plot.data(), t_plot.size());
  mglData t_plotwrited(t_plotwrite.data(), t_plotwrite.size());

  gr.SetRanges(*evals.begin(), *--evals.end(), *t_plot.begin(), *--t_plot.end());
  gr.Plot(evalsd, t_plotd);
  gr.Plot(evalsd, t_plotwrited);
  
  gr.AddLegend("Runtime for plotting", "b");
  gr.AddLegend("Runtime for plotting and writing to eps", "g");
  gr.Legend(1,0);
  gr.Axis();
  gr.Grid("","h");

  gr.Label('x', "number of data points", 0);
  gr.Label('y', "runtime in s", 0);

  gr.WriteEPS("runtimes.eps");
}
  

int main()
{
  measure();
  return 0;
}
