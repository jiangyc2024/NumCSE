template <class Function>
double trapezoidal(Function& f, const double a, const double b, const unsigned N) {
  double I = 0;
  const double h = (b - a)/N; // interval length

  for (unsigned i = 0; i < N; ++i) {
    // rule: T = (b - a)/2 * (f(a) + f(b)),
    // apply on N intervals: [a + i*h, a + (i+1)*h], i=0..(N-1)
    I += h/2*( f(a + i*h) + f(a + (i + 1)*h) );
  }
  return I;
}
