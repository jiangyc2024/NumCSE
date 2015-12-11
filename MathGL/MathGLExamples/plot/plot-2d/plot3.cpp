/*
 * Example on how to do Matlab's plot3(..) in MathGL
 */

# include <Eigen/Dense>
# include <mgl2/mgl.h>
# include <cmath> // atan

void sample(mglGraph* gr)
{
  const double pi = 4*atan(1);
  const std::size_t n = 1000;
  Eigen::VectorXd l = Eigen::VectorXd::LinSpaced(n, 0, 10*pi);
  Eigen::VectorXd x = l.array().sin().matrix();
  Eigen::VectorXd y = l.array().cos().matrix();
  
  Eigen::VectorXd x2 = 0.5*(2*l.array()).sin().matrix();
  Eigen::VectorXd y2 = 0.5*(2*l.array()).cos().matrix();

  mglData xd(x.data(), x.size());
  mglData yd(y.data(), y.size());
  mglData ld(l.data(), l.size());

  mglData x2d(x2.data(), x2.size());
  mglData y2d(y2.data(), y2.size());

  gr->SubPlot(1,1,0,"<_");
  gr->Title("3D Line plot");
  gr->SetRanges(-1,1,-1,1,0,40);
  gr->Rotate(40,60);
  gr->Box();
  gr->Axis();
  gr->Plot(xd,yd,ld);
  gr->Plot(x2d,y2d,ld, "n;");
}

int main()
{
  mglGraph gr;
  sample(&gr);
  gr.WriteEPS("spiral.eps");
  
  return 0;
}
