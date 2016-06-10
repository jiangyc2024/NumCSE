function a = chebexpdir(y)
% Direct computation of coefficients \Blue{$\alpha_j$} in the Chebychev expansion \Blue{$p =
% \sum\limits_{j=0}^{n}\alpha_j T_j$} of \Blue{$p\in\Cp_n$} based on values \Blue{$y_k$},
% \Blue{$k=0,\ldots,n$}, in Chebychev nodes \Blue{$t_k$}, \Blue{$k=0,\ldots,n$}. These values are
% passed in the row vector \texttt{y}.
n = length(y)-1;   % degree of polynomial
t = cos(((0:n)+0.5)/(n+1)*pi); % Chebychev nodes on \Blue{$[-1,1]$}, see \eqref{eq:CHEBNODES}
V = chebpolmult(n,t)'; % Matrix of interpolation conditions, see \eqref{eq:iplse}
% Note: this matrix is in fact orthogonal up to diagonal scaling due to the discrete
% orthogonality of the Chebychev polynomials w.r.t to Chebychev points 
a = (V\reshape(y,n+1,1))'; % Solve linear system (inefficient!)


