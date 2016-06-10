function [alpha, beta] = coeffortho(t,n)
% Vector \texttt{t} passes the points in the definition of the discrete $L^{2}$-inner 
% product, \texttt{n} the maximal index desired
m = numel(t);  % Maximal degree of orthogonal polynomial
alpha(1) = sum(t)/m;
% Initialization of recursion; we store only the values of 
% the polynomials at the points in \Blue{$\Ct$}.
p1 = ones(size(t)); 
p2 = t-alpha(1);
% Main loop
for k=1:min(n-1,m-2)
    p0 = p1;  p1 = p2;
    % 3-term recusion \eqref{eq:onp:rec}, 
    alpha(k+1) = dot(p1,(t.*p1))/norm(p1)^2;
    beta(k) = (norm(p1)/norm(p0))^2;
    p2 = (t-alpha(k+1)).*p1-beta(k)*p0;
end
