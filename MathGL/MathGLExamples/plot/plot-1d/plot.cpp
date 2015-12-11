/*
 * Two short example plots in 1-D
 */

# include <Eigen/Dense>
# include <mgl2/mgl.h>
# include <vector>

void sample_vec(mglGraph* gr)
{
  std::vector<double> v(100);
  for (auto it = v.begin(); it != v.end(); ++it)
    *it = drand48();
  mglData d(v.data(), v.size()); // constructing mglData from vector: works the same for Eigen vectors
  gr->Plot(d);
  gr->Axis();
}

void sample_eig(mglGraph* gr)
{
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(1000, -5, 1);
  Eigen::VectorXd y = (t.array().exp()).matrix(); // manipulate data
  mglData td(t.data(), t.size()); 
  mglData yd(y.data(), y.size());
  gr->SetRanges(t.minCoeff(), t.maxCoeff(), y.minCoeff(), y.maxCoeff()); // set proper ranges for axis
  gr->Plot(td, yd);
  gr->Axis();
  gr->AddLegend("exp(x)", "b");
  gr->Legend();
}

int main()
{
  mglGraph gr_vec, gr_eig;
  sample_vec(&gr_vec);
  gr_vec.WriteEPS("rand.eps");
  sample_eig(&gr_eig);
  gr_eig.WriteEPS("exp.eps");
  return 0;
}
