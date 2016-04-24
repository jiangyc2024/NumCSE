template <class Function>
double simpson(Function& f, const double a, const double b, const unsigned N) {
  double I = 0;
  const double h = (b - a)/N; // intervall length

  for (unsigned i = 0; i < N; ++i) {
    // rule: S = (b - a)/6*( f(a) + 4*f(0.5*(a + b)) + f(b) )
    // apply on [a + i*h, a + (i+1)*h] 
    I += h/6*( f(a + i*h) + 4*f(a + (i+0.5)*h) + f(a + (i+1)*h) );
  }

  return I;
}
