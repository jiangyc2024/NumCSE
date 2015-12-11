# include <natcsi.hpp>
# include <Eigen/Dense>
# include <mgl2/mgl.h>
# include <iostream>
# include <fstream>

int main(){
  // build data, function f(t) = t*exp(sin(t))
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(10, -5, 5);
  Eigen::VectorXd y = (t.array()*t.array().sin().exp()).matrix();
  NatCSI N(t, y);
  const std::size_t n = 200;
  Eigen::VectorXd t_interp = Eigen::VectorXd::LinSpaced(n, t.minCoeff(), t.maxCoeff());
  Eigen::VectorXd y_interp(n);
  for (std::size_t i = 0; i < n; ++i){
    y_interp(i) = N(t_interp(i));
  }
  // prepare data for plotting
  mglData td, yd, td_interp, yd_interp;
  td.Set(t.data(), t.size()); td_interp.Set(t_interp.data(), t_interp.size());
  yd.Set(y.data(), y.size()); yd_interp.Set(y_interp.data(), y_interp.size());

  // plot
  mglGraph gr;
  gr.SetRanges(t_interp.minCoeff(), t_interp.maxCoeff(), y_interp.minCoeff() - 1, y_interp.maxCoeff() + 1);
  gr.Axis();
  gr.Plot(td, yd, " +"); // interpolation data
  gr.Plot(td_interp, yd_interp, "g"); // interpolant
  gr.FPlot("x*exp(sin(x))", "r"); // function
  
  gr.AddLegend("Interpolation data", " +");
  gr.AddLegend("Interpolant", "g");
  gr.AddLegend("Function", "r");
  gr.Legend(1,0);

  gr.WriteEPS("interpolation.eps");


  return 0;
}
