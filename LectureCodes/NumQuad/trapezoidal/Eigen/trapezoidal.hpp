#include <cassert>

// N-interval equidistant trapezoidal rule
template <class Function>
double trapezoidal(Function &f, const double a, const double b,
                   const unsigned N) {
  double I = 0;
  const double h = (b - a) / N; // interval length

  for (unsigned i = 0; i < N; ++i) {
    // rule: T = (b - a)/2 * (f(a) + f(b)),
    // apply on N intervals: [a + i*h, a + (i+1)*h], i=0..(N-1)
    I += h / 2 * (f(a + i * h) + f(a + (i + 1) * h));
  }
  return I;
}

// Alternative implementation of n-point equidistant trapezoidal rule
template <typename Functor>
double equidTrapezoidalRule(Functor &&f, double a, double b, unsigned int n) {
  assert(n>=2);
  const double h = (b - a)/(n-1);
  double t = a + h;
  double s = 0.0;
  for (unsigned int i = 1; i < n-1; t += h, ++i) {
    s += f(t);
  }
  return (0.5 * h * f(a) + 0.5 * h * f(b) + h * s);
}
