/*
 * Short example on legends, fonts and font sizes
 */

# include <mgl2/mgl.h>

void sample(mglGraph* gr)
{
  mglData y0(100), y1a(100), y1b(10), y2(100);
  y0.Modify("sin(pi*x)");
  y1a.Modify("cos(2*pi*x)");
  y1b.Modify("cos(2*pi*x)");
  y2.Modify("exp(-8*x^2)*sin(10*pi*x)");

  gr->AddLegend("sin(\\pi x)", "b");
  gr->SubPlot(2,2,0,"<_"); // title directly over plot
  gr->Title("Default legend & font");
  gr->Plot(y0);
  gr->Box();
  gr->Legend(); // default legend

  gr->AddLegend("cos(\\pi x)", "g");
  gr->SubPlot(2,2,1);
  gr->LoadFont("adventor"); // load font 'adventor'
  gr->Title("Adventor font");
  gr->Plot(y0);
  gr->Plot(y1a);
  gr->Box();
  gr->Legend(1, "Wb#" ); // grey background, blue font and box around it

  gr->AddLegend("exp(-8\\x^2)sin(10\\pi x)", "r");
  gr->SubPlot(2,2,2);
  gr->LoadFont("schola"); // load font 'schola'
  gr->Title("Horizontal & Schola font");
  gr->Plot(y0);
  gr->Plot(y1a);
  gr->Plot(y2);
  gr->Box();
  gr->Legend(0.5, 0, "-#"); // bottom horizontal

  gr->ClearLegend(); // removing all previous legend entries
  gr->SetLegendMarks(2);
  gr->AddLegend("cos(2\\pi x)", "r");
  gr->SubPlot(2,2,3);
  gr->LoadFont("cursor"); // load font 'cursor'
  gr->Title("Cursor font");
  gr->SetRanges(0,1,-1,1); // default is [-1:1],[-1:1] but our data is in [0:1], [-1:1]
  gr->Axis();
  //gr->Plot(y1a,"r");
  gr->Plot(y1b,"rd");
  gr->Legend(0.5,0.5);
}

int main()
{
  mglGraph gr;
  sample(&gr);
  gr.WriteEPS("leg.eps");
  return 0;
}
 
