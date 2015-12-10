/*
 * Comparing polynomial evaluation with
 * horner scheme and straight naive evaluation
 */


# include <iostream>
# include <chrono>
# include <mgl2/mgl.h>
# include <Eigen/Dense>

using std::chrono::high_resolution_clock;
using std::chrono::nanoseconds;
using std::chrono::duration_cast;

using Eigen::VectorXd;

// naive power function
template <typename Scalar>
Scalar pow(Scalar x, int n){
  Scalar res = x;
  for (int i = 1; i < n; ++i)
    res *= x;
  return res;
}

// horner evaluation
template <typename Scalar, typename CoeffVec>
Scalar hornerEval(const CoeffVec& c, const Scalar x){
  Scalar res = c(0);
  for (int i = 1; i < c.size(); ++i)
    res = x*res + c(i);
  return res;
}

// straight evaluation:
// if superbad is true, then a naive power function is used, 
// if superbad is false, the std::pow function is used
template <typename Scalar, typename CoeffVec>
Scalar naiveEval(const CoeffVec& c, const Scalar x, bool superbad = true){
  Scalar res = 0;
  if (superbad){
    for (int i = c.size() - 1; i >= 1; --i)
      res += c(i)*pow(x, i);
  }
  else {
    for (int i = c.size() - 1; i >= 1; --i)
      res += c(i)*std::pow(x, i);
  }
  return res;
}

int main(){
  double x = 0.123;
  const unsigned int N = 100000;
  const unsigned int repeats = 1;

  std::vector<double> evals, horner, naive, supernaive;
  evals.reserve(N);
  horner.reserve(N);
  naive.reserve(N);
  supernaive.reserve(N);

  // testing for: horner scheme - naive with naive power-function - naive with efficient power function
  for (unsigned int d = 2; d < N; d *=2){
    // random coefficients for polynom
    Eigen::VectorXd c = Eigen::VectorXd::Random(d);
    // saving degrees for which we evaluated
    evals.push_back(d);
    double buffer;

    // horner scheme
    auto ht = high_resolution_clock::now();
    for (unsigned int i = 0; i < repeats; ++i)
      buffer = hornerEval(c, x);
    horner.push_back(duration_cast<nanoseconds>(high_resolution_clock::now() - ht).count()/double(1e9)); // normalize data to seconds
    std::cout << buffer; // print to make sure the function doesnt get optimized away
    
    // naive evaluation with naive power function     
    auto nt = high_resolution_clock::now();
    for (unsigned int i = 0; i < repeats; ++i)
      buffer = naiveEval(c, x);
    supernaive.push_back(duration_cast<nanoseconds>(high_resolution_clock::now() - nt).count()/double(1e9));
    std::cout << buffer;
    
    // naive evaluation with std::pow
    auto nt2 = high_resolution_clock::now();
    for (unsigned int i = 0; i < repeats; ++i)
      buffer = naiveEval(c, x, false);
    naive.push_back(duration_cast<nanoseconds>(high_resolution_clock::now() - nt2).count()/double(1e9));
    std::cout << buffer;
  }

  // preparing data for plot
  mglData evalsd(evals.data(), evals.size());
  mglData hornerd(horner.data(), horner.size());
  mglData supernaived(supernaive.data(), supernaive.size());
  mglData naived(naive.data(), naive.size());

  // plotting results
  mglGraph gr;
  gr.SubPlot(1,1,0,"<_"); // with this trick the title will appear directly over the plot
  gr.SetFontSizePT(6);
  gr.Title("Polynomial evaluation"); // only problem: title must be shorter/smaller
  gr.SetRanges(evalsd.Minimal(), evalsd.Maximal(), hornerd.Minimal(), naived.Maximal());
  gr.SetFunc("lg(x)","lg(y)");
  
  gr.Axis();
  gr.Grid("","h");
 
  // plot data
  gr.Plot(evalsd, hornerd, "g#+");
  gr.Plot(evalsd, supernaived, "r+");
  gr.Plot(evalsd, naived, "b^");

  // using FPlot for comparison-lines
  gr.FPlot("x/3e7","k:");
  gr.FPlot("x^2/2e8","k;");

  gr.AddLegend("Horner", "g#+");
  gr.AddLegend("Naive w/ naive power", "r+");
  gr.AddLegend("Naive w/ std::pow", "b^");
  gr.AddLegend("O(n^2)", "k;");
  gr.AddLegend("O(n)", "k:");
  gr.Legend(0,1);

  gr.WriteEPS("runtimes.eps");
  
  return 0;
}
  
  
