# include <interpol.hpp>
# include <algorithm>
# include <numeric> // accumulate
# include <functional> // mulitplies

Interpol::Interpol(Eigen::VectorXd& t, Eigen::VectorXd& y)
{
  assert(t.size() == y.size()); // should be implemented with try and catch
  const long n = t.size() - 1; // t = t_0,... ,t_n -> size = n + 1
  Eigen::VectorXd lambda = Eigen::VectorXd::Zero(n + 1);

  for (long i = 0; i < n + 1; ++i){
    std::vector<double> T;
    T.reserve(n);
    for (long j = 0; j < n + 1; ++j){
      if (i != j)
        T.push_back(t(i) - t(j)); 
    }

    // following line computes the product of the elements in T
    double T_prod = std::accumulate(T.begin(), T.end(), 1.0, std::multiplies<double>());
    lambda(i) = 1./T_prod;
  }
  t_ = t;
  y_ = y;
  lambda_ = lambda;
}

double Interpol::operator()(double x){
  auto ind_ptr = std::find(t_.data(), t_.data() + t_.size(), x);
  int ind = ind_ptr - t_.data(); // get index as number
  // if x has the same value as a node we must avoid division by 
  // zero and return the value at the node
  if (ind_ptr != t_.data() + t_.size())
    return y_(ind);
  // else use baryzentric formula
  Eigen::VectorXd mu = (lambda_.array()/(x - t_.array())).matrix();
  double result = (mu.array()*y_.array()).sum()/mu.sum();
  return result;
}
