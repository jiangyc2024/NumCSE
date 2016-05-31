class ChebInterp {
    private: 
        // various internal data describing Chebychev interpolating polynomial \Blue{$p$}
    public: 
       // Constructor taking function \Blue{$f$} and degree \Blue{$n$} as arguments
       template <typename Function>
       ChebInterp(const Function &f, unsigned int n);
       // Evaluation operator: \Blue{$y_{j} = p(x_{j})$}, \Blue{$j=1,\ldots,m$} (\Blue{$m$} ``large'')
       double eval(const vector<double> &x, vector<double> &y) const;
};
