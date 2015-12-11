/*
 * 2-D Plot examples in MathGL
 */

# include <iostream>
# include <Eigen/Dense>
# include <mgl2/mgl.h>
# include "grid.hpp"


using namespace Eigen;

void sample(mglGraph* gr1, mglGraph* gr2, mglGraph* gr3, mglGraph* gr4)
{
  VectorXd l = VectorXd::LinSpaced(100,-4,4);
  std::pair<MatrixXd, MatrixXd> Grid = meshgrid(l,l); // "selfmade" meshgrid function in grid.hpp
  MatrixXd Z = (2*(Grid.first.array()).sin()*(2.5*Grid.second.array()).cos()).matrix();
  mglData Xd(Grid.first.rows(), Grid.first.cols(), Grid.first.data());
  mglData Yd(Grid.second.rows(), Grid.second.cols(), Grid.second.data());
  mglData Zd(Z.rows(), Z.cols(), Z.data());
 
  // meshplot example
  gr1->Title("Mesh plot");
  gr1->SetRanges(-1,1,-1,1,-2,2);
  gr1->Rotate(50,60);
  gr1->Box();
  gr1->Mesh(Xd,Yd,Zd);
  gr1->Colorbar(); // Default colorbar is right outside the box
  gr1->Axis();

  // surf example
  gr2->Title("Surf plot");
  gr2->SetRanges(-4,4,-4,4,-2,2);
  gr2->Rotate(50,60);
  gr2->Box();
  gr2->Grid("xyz","h");
  // gr2->Surf(Zd) also gives *some* results but to have control its better to also give Xd and Yd as arguments
  gr2->Surf(Xd,Yd,Zd); 
  gr2->Axis();

  // bad example: Ranges are not automatically fitted to the data
  gr3->Title("Bad auto ranges");
  gr3->Rotate(50,80);
  gr3->Box();
  gr3->Mesh(Zd);
  gr3->Axis();

  // contour example
  gr4->Title("Contour plot");
  gr4->SetRange('c',-2,2); // 'c' for colobar
  gr4->Cont(Zd); // ContF for filled contour plot
  gr4->Colorbar("I<"); // < = right, I = next to the box
}

int main()
{
  mglGraph gr1, gr2, gr3, gr4;
  sample(&gr1, &gr2, &gr3, &gr4);
  gr1.WriteEPS("mesh.eps");
  gr2.WriteEPS("surf.eps");
  gr3.WriteEPS("bad.eps");
  gr4.WriteEPS("cont.eps");
  return 0;
}
