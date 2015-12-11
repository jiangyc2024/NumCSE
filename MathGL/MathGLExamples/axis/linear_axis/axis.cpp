/*
 * Short example on axis, scales and labels
 */

# include <mgl2/mgl.h>

void sample(mglGraph* gr)
{
  mglData y0(200), y1(200), y2(200); // gives data in [0,1]
  y0.Modify("exp(-8*x^2)*sin(10*pi*x) + 1");
  y1.Modify("exp(-8*x^2)*sin(5*pi*x) + 1");
  y2.Modify("exp(-8*x^2)*cos(5*pi*x) + 1");
  
  // ------------------ first plot ------------------- //
  gr->SubPlot(2,2,0);
  gr->Title("Default axis");
  gr->Axis();
  gr->Puts(0,-0.5,"Scales are not necessarily sensible in default axis!", "#");
  gr->Plot(y0);
  gr->Plot(y1);
  gr->Plot(y2);

  // ----------------- second plot ------------------- //
  gr->SubPlot(2,2,1);
  gr->Title("Adjusted axis");
  gr->SetRange('y',0,2);
  gr->SetOrigin(NAN,NAN); // with this options MathGL tries to find best origin
  gr->Axis();
  gr->Label('x', "x-Axis", 1); // alignment: -1 = left, 0 = center, 1 = right
  gr->Label('y', "y-Axis", 0);
  gr->Plot(y0);
  gr->Plot(y1);
  gr->Plot(y2);

  // ------------------ third plot ------------------- //
  gr->SubPlot(2,2,2);
  gr->Title("Multiple scales");
  // setting first axis
  gr->SetRange('x',-2,2); // these two lines are equivalent to:
  gr->SetRange('y',-1,1); // gr->SetRanges(-2,2,-1,1);
  gr->SetOrigin(0,0);
  gr->Axis();
  // plotting in scale defined by first axis
  gr->FPlot("x^4-x^2");
  
  // setting second axis
  gr->SetRanges(-2,2,0,2);
  gr->SetOrigin(-2,0);
  gr->Axis("y","r"); // "y" -> plotting only y-axis, "r" -> in red
  // plotting in scale defined by second axis
  gr->FPlot("x^2","r");

  // adding legend
  gr->AddLegend("\\x^4 - \\x^2", "b");
  gr->AddLegend("\\x^2", "r");
  gr->Legend(1,0,"#");

  gr->SubPlot(2,2,3);
  gr->SetFontSizePT(2);
  gr->Title("_{Tiny}, @{Small}, Normal, \\big{Big}");
  gr->Puts(mglPoint(0,0.7), "@{Sample}");
  gr->Puts(mglPoint(0,0.5), "Sample");
  gr->Puts(mglPoint(0,0.2), "\\big{Sample}");

  return;
}

int main()
{
  mglGraph gr;
  sample(&gr);
  gr.WriteEPS("scales.eps");
  return 0;
}
