function [a,b] = trigpolycoeff(t,y)
% Computes expansion coefficients of trigonometric polyonomials \eqref{eq:trigpreal}
% \texttt{t}: row vector of nodes \Blue{$t_{0},\ldots,t_n\in[0,1[$}
% \texttt{y}: row vector of data \Blue{$y_{0},\ldots,y_{n}$}
% return values are column vectors of expansion coefficients \Blue{$\alpha_j$}, \Blue{$\beta_j$}
N = length(y); if (mod(N,2)~=1), error('#pts odd!'); end
n = (N-1)/2;
M = [cos(2*pi*t'*(0:n)), sin(2*pi*t'*(1:n))];
c = M\y';
a = c(1:n+1); b = c(n+2:end);
 
