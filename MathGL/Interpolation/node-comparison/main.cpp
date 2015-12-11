# include <iostream>
# include <cmath>
# include <vector>
# include "interpol.hpp" // interpol already includes Eigen/Dense
# include <mgl2/mgl.h>

double interpol(Eigen::VectorXd t, mglGraph* gr = 0)
{
  Eigen::VectorXd y = (1/(1 + t.array()*t.array())).matrix();

  Interpol intp(t, y);

  long N = 500;
  Eigen::VectorXd t_intp = Eigen::VectorXd::LinSpaced(N, -5, 5);
  Eigen::VectorXd y_intp(t_intp.size());
  for (long i = 0; i < N; ++i)
    y_intp(i) = intp(t_intp(i));

  if (gr != 0){
    // plot interpolation
    mglData td(t.data(), t.size());
    mglData yd(y.data(), y.size());
    mglData td_intp(t_intp.data(), t_intp.size());
    mglData yd_intp(y_intp.data(), y_intp.size());

    gr->SetRanges(-5.1,5.1,-0.5,1.5);
    gr->Axis();
    gr->Grid("","h");
    gr->FPlot("1/(1+x^2)","r");
    gr->Plot(td, yd, " no");
    gr->Plot(td_intp, yd_intp, "b;");

    gr->AddLegend("\\(1 + x^2)^{-1}", "r");
    gr->AddLegend("Interpolant", "b;");
    gr->AddLegend("Interpolation data", " no");
    gr->Legend();
  }
  // compute maximal error
  double max_err = ((1/(1 + t_intp.array()*t_intp.array())).matrix() - y_intp).maxCoeff();
  return max_err;
}

int main()
{

  /** Plotting for n = 10 nodes **/
  long n = 10;
  // equidistant nodes
  Eigen::VectorXd t_equi = Eigen::VectorXd::LinSpaced(n + 1, -5, 5);
  // chebychev nodes
  Eigen::VectorXd t_cheb(n + 1);
  const double pi = 4*atan(1);
  for (int i = 0; i <= n; ++i)
    t_cheb(i) = -5 + 5*(cos((2*i + 1)*pi/(2*n + 2)) + 1);

  mglGraph gr_equi, gr_cheb;
  gr_equi.Title("Equidistant nodes");
  interpol(t_equi, &gr_equi);
  gr_equi.WriteEPS("intp-equi.eps");

  gr_cheb.Title("Chebychev nodes");
  interpol(t_cheb, &gr_cheb);
  gr_cheb.WriteEPS("intp-cheb.eps");

  /** Plotting error **/
  long N = 100;

  
  // computing error
  std::vector<double> err_equi, err_cheb, evals;
  err_equi.reserve(N/2);
  err_cheb.reserve(N/2);
  evals.reserve(N/2);
  
  for(long n = 10; n < N; n += 2){
    evals.push_back(n);
    // equidistant nodes
    Eigen::VectorXd t_equi = Eigen::VectorXd::LinSpaced(n + 1, -5, 5);
    
    // chebychev nodes
    Eigen::VectorXd t_cheb(n + 1);
    for (int i = 0; i <= n; ++i)
      t_cheb(i) = -5 + 5*(cos((2*i + 1)*pi/(2*n + 2)) + 1);

    err_equi.push_back(interpol(t_equi));
    err_cheb.push_back(interpol(t_cheb));
  }

  // preparing data for plot
  mglData evals_data(evals.data(), evals.size());
  mglData cheb_data(err_cheb.data(), err_cheb.size());
  mglData equi_data(err_equi.data(), err_equi.size());

  // plotting
  mglGraph gr_error;
  gr_error.SubPlot(1,1,0,"<_");
  gr_error.Title("Error");
  gr_error.SetRanges(10, N, 1e-8, 1e2);
  gr_error.SetFunc("","lg(y)");
  gr_error.Label('x', "x", 0);
  gr_error.Label('y', "log(y)", 0);

  gr_error.Plot(evals_data, equi_data, "r-+");
  gr_error.Plot(evals_data, cheb_data, "b-+");

  gr_error.AddLegend("Equidistant nodes (L^{\\infty}-norm)", "r-+");
  gr_error.AddLegend("Chebychev nodes (L^{\\infty}-norm)", "b-+");
  gr_error.Legend(1,0.8);

  gr_error.Axis();
  gr_error.Grid("","h");

  gr_error.WriteEPS("error.eps");
    

  return 0;
}
