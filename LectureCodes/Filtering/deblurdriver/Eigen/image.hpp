# include <Eigen/Dense>
# include <mgl2/mgl.h>

void image(Eigen::MatrixXd A, const std::string& title, const std::string& epsfile) {
  A /= A.maxCoeff(); // normalize data
  mglGraph gr;
  // this line removes too large white margins of the graphic
  // and leaves space for a title 
  gr.SubPlot(1,1,0,"^"); 
  gr.Title(title.c_str()); // write title
  // intitialize MathGL data with Eigen Matrix
  mglData Ad(A.cols(), A.rows(), A.data());

  // set range of colorbar ('c')
  gr.SetRange('c', 0, 1);
  gr.Box(); // draw box around graph
  gr.SetFontSizePT(6);
  gr.Colorbar("kw"); // plot colorbar
  gr.Tile(Ad, "kw"); // plot tiles
  gr.WriteEPS(epsfile.c_str()); // save plot
}
