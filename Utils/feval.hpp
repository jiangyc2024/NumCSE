# ifndef feval_hpp
# define feval_hpp

# include <Eigen/Dense>

using vector_t = Eigen::VectorXd;

// feval evaluates the scalar function f at all elements of the vector x
template <class Function>
void feval(const Function& f, const vector_t& x, vector_t& y) {
  y = vector_t(x.size());
  for (int i = 0; i < x.size(); ++i) {
    y(i) = f(x(i));
  }
}

// same as feval above but explicitly returns f(x)
template <class Function>
vector_t feval(const Function& f, const vector_t& x) {
  vector_t y;
  feval(f, x, y);
  return y;
}

# endif
