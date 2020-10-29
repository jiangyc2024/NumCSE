class PolyEval {
private:
  // evaluation point and various internal data describing the polynomials
public:
  // Constructor taking the evaluation point as argument
  PolyEval(double x);
  // Add another data point and update internal information
  void addPoint(double t, double y);
  // Value of current interpolating polynomial at \texttt{x}
  double eval(void) const;
};
