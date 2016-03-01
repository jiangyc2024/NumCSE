function I = quadformula(f,c,w)
% Generic numerical quadrature routine implementing \eqref{eq:quadform}: 
% \texttt{f} is a handle to a function of type \texttt{f = @(x) ....}
% \texttt{c}, \texttt{w} pass quadrature nodes \Blue{$\qn_j\in [a,b]$}, and weights \Blue{$\qw_j>0$}
n = length(c); If (length(w) ~= n), error('#weights != #nodes'); end
I = 0; for j=1:n, I = I + w(j)*f(c(j)); end
  
