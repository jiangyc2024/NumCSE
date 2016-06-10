function a = chebexp(y)
% Efficiently compute coefficients \Blue{$\alpha_j$} in the Chebychev expansion 
% \Blue{$p = \sum\limits_{j=0}^{n}\alpha_j T_j$} of \Blue{$p\in\Cp_n$} based on values \Blue{$y_k$},
% \Blue{$k=0,\ldots,n$}, in Chebychev nodes \Blue{$t_k$}, \Blue{$k=0,\ldots,n$}. These values are
% passed in the row vector \texttt{y}.
n = length(y)-1;   % degree of polynomial
% create vector \Blue{$\Vz$} by wrapping and componentwise scaling
z = exp(-pi*i*n/(n+1)*(0:2*n+1)).*[y,y(end:-1:1)]; % r.h.s. vector 
c = ifft(z);  % Solve linear system \eqref{eq:tpiplse} with effort \Blue{$O(n\log n)$}
b = real(exp(0.5*pi*i/(n+1)*(-n:n+1)).*c); % recover \Blue{$\beta_j$}, see \eqref{eq:tpiplse}
a = [b(n+1),2*b(n+2:2*n+1)]; % recover \Blue{$\alpha_j$}, see \eqref{eq:tpdft}



