///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): 
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#define NDEBUG true //NOLINT(cppcoreguidelines-macro-usage)
#include <Eigen/Dense>
#include <chrono>
#include <cstdint>
#include <iostream>

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;
using index_t = Eigen::Index;
using value_t = Eigen::MatrixXd::Scalar;

/** Timing of row and column oriented access for Eigen */
/* SAM_LISTING_BEGIN_1 */
void rowcolaccesstiming()
{
  constexpr size_t K = 3; // Number of repetitions  
  constexpr index_t N_min = 5; // Smallest matrix size 32
  constexpr index_t N_max = 13; // Scan until matrix size of 8192
  index_t n = (1UL << static_cast<size_t>(N_min));
  Eigen::MatrixXd times(N_max-N_min+1,3);
  
  for(index_t l=N_min; l<= N_max; l++, n*=2) {
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
    value_t t1 = 1000.0;
    for(size_t k=0;k<K;k++) {
      auto tic =  high_resolution_clock::now();
      for(index_t j=0; j < n-1; j++) {
        A.row(j+1) -= A.row(j); // row access
      }
      auto toc = high_resolution_clock::now();
      const value_t t = static_cast<value_t>(duration_cast<microseconds>(toc-tic).count())/1E6;
      t1 = std::min(t1,t);
    }
    value_t t2 = 1000.0;
    for(size_t k=0;k<K;k++) {
      auto tic =  high_resolution_clock::now();
      for(index_t j=0; j < n-1; j++) {
        A.col(j+1) -= A.col(j); //column access
      }
      auto toc = high_resolution_clock::now();
      const value_t t = static_cast<value_t>(duration_cast<microseconds>(toc-tic).count())/1E6;
      t2 = std::min(t2,t);
    }
    times(l-N_min,0) = static_cast<value_t>(n);
    times(l-N_min,1) = t1;
    times(l-N_min,2) = t2;
  }
  std::cout << times << std::endl;
}
/* SAM_LISTING_END_1 */

// Overhead of returning by value
// First version return by reference
void retMatRef(const Eigen::VectorXd &v,Eigen::MatrixXd &R)
{
  using index_t = typename Eigen::VectorXd::Index;
  const index_t n = v.size();
  R = (v*v.transpose() + Eigen::MatrixXd::Identity(n,n));
}

Eigen::MatrixXd retMatVal(const Eigen::VectorXd &v)
{
  using index_t = typename Eigen::VectorXd::Index;
  const index_t n = v.size();
  return (v*v.transpose() + Eigen::MatrixXd::Identity(n,n));
}

void retmattiming()
{
  std::cout << "Timing return by reference vs. return by value" << std::endl;
  constexpr size_t K = 3; // Number of repetitions
  constexpr index_t N_min = 4; // Smallest matrix size 16
  constexpr index_t N_max = 13; // Scan until matrix size of 8192
  Eigen::MatrixXd res(N_max-N_min+1,4);
  index_t n = (1UL << static_cast<size_t>(N_min));
  
  for(index_t l=N_min; l<= N_max; l++, n*=2) {
    Eigen::VectorXd v = Eigen::VectorXd::LinSpaced(n,0.0,1.0);
    Eigen::MatrixXd M(n,n);
    value_t t1 = 10000.0;
    for(size_t k=0;k<K;k++) {
      auto tic =  high_resolution_clock::now();
      retMatRef(v,M); M(n-1,n-1) += static_cast<value_t>(k);
      auto toc = high_resolution_clock::now();
      const value_t t = static_cast<value_t>(duration_cast<microseconds>(toc-tic).count())/1E6;
      t1 = std::min(t1,t);
      v *= 1.5;
    }
    value_t t2 = 10000.0;
    for(size_t k=0;k<K;k++) {
      auto tic =  high_resolution_clock::now();
      M = retMatVal(v);  M(n-1,n-1) += static_cast<value_t>(k); //NOLINT(clang-analyzer-core.uninitialized.Assign) //bug in NDEBUG version of Eigen
      auto toc = high_resolution_clock::now();
      const value_t t = static_cast<value_t>(duration_cast<microseconds>(toc-tic).count())/1E6;
      t2 = std::min(t2,t);
      v *= 1.5;
    }
    std::cout << "n = " << n << ", t(Ref) = " << t1 << ", t(Val) " << t2 << std::endl;
    res(l-N_min,0) = static_cast<value_t>(n);
    res(l-N_min,1) = t1;
    res(l-N_min,2) = t2;
    res(l-N_min,3) = t1/t2; 
  }
  std::cout << "n t(Ref) t(val) tr/tv" << std::endl << res << std::endl;
}

int main(int argc,char **argv) {
  int code = 0;
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <selection>" << std::endl; //NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    code = -1;
  }
  else {
    const int64_t sel = std::strtol(argv[1], nullptr, 10); //NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    switch (sel) {
    case 1: { rowcolaccesstiming(); break; }
    case 2: { retmattiming(); break; }
    default: { std::cerr << "Invalid selection" << std::endl; code = -1; }
    }
  }
  return code;
}
