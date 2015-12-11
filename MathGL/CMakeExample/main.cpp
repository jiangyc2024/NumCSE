# include <mgl2/mgl.h>
# include <mylibrary.hpp>

int main()
{
  std::vector<double> v = sample();
  mglData d;
  d.Set(v);
  mglGraph gr;
  gr.SetRanges(0,10,0,100);
  gr.Axis();
  gr.Plot(d, "r-+");
  gr.WriteEPS("testplot.eps");
  return 0;
}
