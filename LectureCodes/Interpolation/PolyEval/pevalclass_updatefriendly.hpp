class PolyEval {
private:
  std::vector<double> t;  // Interpolation nodes
  std::vector<double> y;  // Coefficients in Newton representation
public:
  PolyEval(); // Idle constructor
  void addPoint(double t, double y); // Add another data point
  // evaluate value of current interpolating polynomial at \Blue{$x$}, 
  double operator() (double x) const;
private:
  // Update internal representation, called by \texttt{addPoint()}
  void divdiff();
};

