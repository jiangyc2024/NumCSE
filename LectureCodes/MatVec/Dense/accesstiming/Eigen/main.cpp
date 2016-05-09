#define NDEBUG true
#include <Eigen/Dense>
#include <iostream>
#include <chrono>
using namespace std::chrono;

/** Timing of row and column oriented access for Eigen */
void rowcolaccesstiming(void)
{
  const int K = 3; // Number of repetitions
  const int N_min = 4; // Smalles matrix size 16
  const int N_max = 13; // Scan until matrix size of 8192
  unsigned long n = (1L << N_min); 
  
  for(int l=4; l<= N_max; l++, n*=2) {
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
    double t1 = 1000.0;
    for(int k=0;k<K;k++) {
      auto tic =  high_resolution_clock::now();
      for(int j=0; j < n-1; j++) A.row(j+1) -= A.row(j); // row access
      auto toc = high_resolution_clock::now();
      double t = (double)duration_cast<microseconds>(toc-tic).count()/1E6;
      t1 = std::min(t1,t);
    }
    double t2 = 1000.0;
    for(int k=0;k<K;k++) {
      auto tic =  high_resolution_clock::now();
      for(int j=0; j < n-1; j++) A.col(j+1) -= A.col(j); //column access
      auto toc = high_resolution_clock::now();
      double t = (double)duration_cast<microseconds>(toc-tic).count()/1E6;
      t2 = std::min(t2,t);
    }
    std::cout << "n = " << n << ", t(row) = " << t1 << ", t(col) " << t2 << std::endl;
  }
}

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

void retmattiming(void)
{
  std::cout << "Timing return by reference vs. return by value" << std::endl;
  const int K = 3; // Number of repetitions
  const int N_min = 4; // Smalles matrix size 16
  const int N_max = 13; // Scan until matrix size of 8192
  Eigen::MatrixXd res(N_max-N_min+1,4);
  unsigned long n = (1L << N_min); 
  
  for(int l=N_min; l<= N_max; l++, n*=2) {
    Eigen::VectorXd v = Eigen::VectorXd::LinSpaced(n,0.0,1.0);
    Eigen::MatrixXd M(n,n);
    double t1 = 10000.0;
    for(int k=0;k<K;k++) {
      auto tic =  high_resolution_clock::now();
      retMatRef(v,M); M(n-1,n-1) += k;
      auto toc = high_resolution_clock::now();
      double t = (double)duration_cast<microseconds>(toc-tic).count()/1E6;
      t1 = std::min(t1,t);
      v *= 1.5;
    }
    double t2 = 10000.0;
    for(int k=0;k<K;k++) {
      auto tic =  high_resolution_clock::now();
      M = retMatVal(v);  M(n-1,n-1) += k;
      auto toc = high_resolution_clock::now();
      double t = (double)duration_cast<microseconds>(toc-tic).count()/1E6;
      t2 = std::min(t2,t);
      v *= 1.5;
    }
    std::cout << "n = " << n << ", t(Ref) = " << t1 << ", t(Val) " << t2 << std::endl;
    res(l-N_min,0) = n; res(l-N_min,1) = t1; res(l-N_min,2) = t2;  res(l-N_min,3) = t1/t2; 
  }
  std::cout << "n t(Ref) t(val) tr/tv" << std::endl << res << std::endl;
}

int main(int argc,char **argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <selection>" << std::endl;
    return(-1L);
  }
  else {
    const int sel = atoi(argv[1]);
    switch (sel) {
    case 1: { rowcolaccesstiming(); break; }
    case 2: { retmattiming(); break; }
    default: { std::cerr << "Invalid selection" << std::endl; exit(-1L); }
    }
  }
  return 0;
}
