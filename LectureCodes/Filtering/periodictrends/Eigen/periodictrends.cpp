# include <cmath> // floor
# include <vector>
# include <string> // string, stod
# include <fstream> // getline
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
# include <figure/figure.hpp>
using Eigen::VectorXd;
using Eigen::VectorXcd;

VectorXd read(const std::string file) {
  std::fstream fs(file); // filestream
  std::vector<double> data; // making use of push_back
  std::string buffer; // saving content of current line
  while (std::getline(fs, buffer)) { // 'for every line'
    data.push_back(std::stod(buffer)); // save number
  }
  // finally convert to eigen vector
  Eigen::Map<VectorXd> x(data.data(), data.size()); 
  return x;
}

int main() {
/* SAM_LISTING_BEGIN_0 */
  VectorXd x = read("trend.dat");
  const int n = x.size();
  
  mgl::Figure trend;
  trend.plot(x, "r");
  trend.title("Google: 'Vorlesungsverzeichnis'");
  trend.grid();
  trend.xlabel("week (1.1.2004-31.12.2010)");
  trend.ylabel("relative no. of searches");
  trend.save("searchdata");

  Eigen::FFT<double> fft;
  VectorXcd c = fft.fwd(x);
  VectorXd p = c.cwiseProduct(c.conjugate()).real()
                .segment(2, std::floor((n+1.)/2));

  // plot power spectrum
  mgl::Figure fig; fig.plot(p,"m"); fig.grid();
  fig.title("Fourier spectrum");
  fig.xlabel("index j of Fourier component");
  fig.ylabel("|c_j|^2");

  // mark 4 highest peaks with red star
  VectorXd p_sorted = p; 
  // sort descending
  std::sort(p_sorted.data(), p_sorted.data() + p.size(), std::greater<double>());
  for (int i = 0; i < 4; ++i) {
    // get position of p_sorted[i]
    int idx = std::find(p.data(), p.data()+p.size(), p_sorted[i]) - p.data();
    std::vector<double> ind = { double(idx + 1) }; // save in vectors
    std::vector<double> val = { p_sorted[i] }; // to be able to plot
    fig.plot(ind, val, " r*");
  }
  fig.save("fourierdata");
/* SAM_LISTING_END_0 */

  return 0;
}
