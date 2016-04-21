class Interpolant {
  private:
    // Various internal data describing \Blue{$f$}
    // Can be the coefficients of a basis representation \eqref{eq:funcrep}
  public:
    // Constructor: computation of coefficients \Blue{$c_j$} of representation \eqref{eq:funcrep}
    @\com{Interpolant}@(const vector<double>& t, const vector<double>& y);
    // Evaluation operator for interpolant \Blue{$f$}
    double @\com{operator()} (double t) const;
};
