class PolyEval {
  private:
    // evaluation point and various internal data describing the polynomials
  public:
    // Idle Constructor 
    PolyEval();
    // Add another data point and update internal information
    void addPoint(t,y); 
    // Evaluation of current interpolating polynomial at \texttt{x}
    Eigen::VectorXd operator () (const Eigen::VectorXd &x) const; 
};  
