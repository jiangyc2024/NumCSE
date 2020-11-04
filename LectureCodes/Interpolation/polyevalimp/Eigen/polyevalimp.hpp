PolyEval::PolyEval() {}

void PolyEval::addPoint(double td, double yd) {
  t.push_back(td);
  y.push_back(yd);
  int n = t.size();
  for (int j = 0; j < n - 1; j++)
    y[n - 1] = ((y[n - 1] - y[j]) / (t[n - 1] - t[j]));
}

double PolyEval::operator()(double x) const {
  double s = y.back();
  for (int i = y.size() - 2; i >= 0; --i)
    s = s * (x - t[i]) + y[i];
  return s;
}
