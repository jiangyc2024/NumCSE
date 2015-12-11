/*
 * Short example on logarithmic scales
 */

# include <mgl2/mgl.h>

void sample(mglGraph* gr)
{
  gr->SubPlot(2,2,0);
  gr->Title("Semilog w/ grid");
  // setting semilog scale and plotting
  gr->SetRanges(0,2,0.01,10);
  gr->SetFunc("","lg(y)"); // set semilog-y
  gr->FPlot("exp(x) - 1");
  gr->Grid(); // normal grid
  gr->Axis();
  gr->Label('y',"log(y)",0);
  gr->Label('x',"x", 0);

  // adding legend
  gr->AddLegend("exp(x)", "b");
  gr->Legend(1,0);
  gr->ClearLegend();

  gr->SubPlot(2,2,1);
  gr->Title("Loglog");
  // setting loglog scale and plotting
  gr->SetRanges(0.1,1,0.0001,1);
  gr->SetFunc("lg(x)","lg(y)"); // set loglog plot
  gr->FPlot("x^1.6", "b:"); // linestyle: dashed = .......
  gr->FPlot("x^2", "gi"); // linestyle: dash-dotted = _._._. 
  gr->FPlot("x^3", "r;"); // linestyle: dashed = ------
  gr->Axis();
  gr->Label('x', "log(x)", 0);
  gr->Label('y', "log(y)", 0);

  // adding legend
  gr->AddLegend("\\x^{1.6}", "b:");
  gr->AddLegend("\\x^2", "gi");
  gr->AddLegend("\\x^3", "r;");
  gr->Legend(0,1,"#");
  gr->ClearLegend();

  gr->SubPlot(2,2,2);
  gr->Title("Minor grid");
  gr->SetFunc("",""); // unset loglog plot
  gr->SetRanges(-1,1,-5,5);
  gr->FPlot("1/(x^2 + 1)");
  gr->Axis();
  gr->Grid("!","h="); // ! for minor grid, h for gray, = for dashed
  
  // adding legend 
  gr->AddLegend("\\1/(x^2 + 1)", "b");
  gr->Legend(0.5,0);
  gr->ClearLegend();

  gr->SubPlot(2,2,3);
  gr->Title("Minor grid in loglog");
  gr->SetRanges(0.1,1,0.1,1);
  gr->SetFunc("lg(x)","lg(y)");
  gr->FPlot("x*exp(x)","q");
  gr->Grid("!","h="); // ! for minor grid, h for gray, = for dashed
  gr->Axis();
  
  // adding legend
  gr->AddLegend("xexp(x)","q");
  gr->Legend(1,0);
}

int main()
{
  mglGraph gr;
  sample(&gr);
  gr.WriteEPS("logscales.eps");
  return 0;
}
