function x = lsqtotal(A,b);
% computes only solution \Blue{$\Vx$} of fitted consistent LSE
[m,n]=size(A);
[U, Sigma, V] = svd([A,b]); % see \eqref{tlsq:1}
s = V(n+1,n+1);
if s == 0, error('No solution'); end
x = -V(1:n,n+1)/s; % see \eqref{tlsq:3}