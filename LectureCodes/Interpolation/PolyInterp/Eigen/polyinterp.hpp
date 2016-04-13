class PolyInterp {
  private:
    // various internal data describing \Blue{$p$}
    Eigen::VectorXd t;
  public:
    // Constructors taking node vector \Blue{$(t_{0},\ldots,t_{n})$} as argument
  PolyInterp(const Eigen::VectorXd &_t);
  template <typename SeqContainer> 
           PolyInterp(const SeqContainer &v); 
  // Evaluation operator for data \Blue{$(y_{0},\ldots,y_{n})$}; computes 
  // \Blue{$p(x_{k})$} for \Blue{$x_{k}$}'s passed in \texttt{x}
  Eigen::VectorXd eval(const Eigen::VectorXd &y,
                       const Eigen::VectorXd &x) const;
 };  
