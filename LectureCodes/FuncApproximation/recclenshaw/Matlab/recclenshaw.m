function y = recclenshaw(a,x)
% Recursive evaluation of a polynomial \Blue{$p = \sum_{j=1}^{n+1}a_j T_{j-1}$} at point \texttt{x}, see \eqref{eq:cstr}
% The coefficients \Blue{$a_j$} have to be passed in a row vector.
n = length(a)-1;
if (n<2), y = a(1)+x*a(2);
else 
  y = recclenshaw([a(1:n-2), a(n-1)-a(n+1), a(n)+2*x*a(n+1)],x);
end
