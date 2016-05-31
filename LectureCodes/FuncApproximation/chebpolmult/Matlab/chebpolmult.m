function V = chebpolmult(d,x)
% Computes the values of the Chebychev polynomials \Blue{$T_{0},\ldots,T_{d}$}, \Blue{$d\geq 2$}
% at points passed in \texttt{x} using the 3-term recursion \eqref{eq:ChebychevRecursion}.
% The values \Blue{$T_k(x_j)$}, are returned in row \Blue{$k+1$} of \texttt{V}.
V = ones(1,length(x)); % \Blue{$T_0 \equiv 1$}
V = [V; x];            % \Blue{$T_1(x) = x$}
for k=1:d-1
  V = [V; 2*x.*V(end,:)-V(end-1,:)]; % \Magenta{3-term recursion}
end
