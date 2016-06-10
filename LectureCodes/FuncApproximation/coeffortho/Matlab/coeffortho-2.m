function [alpha, beta, b, P] = coeffortho(t,y,n)
% COEFFORTHO computes the coefficients of the 3-term 
% recurrence for the orthogonal polynomials to the 
% set of points t_i. P is the orthogonal matrix 
% which contains as columns  the values p_i(t_j) of the 
% orthogonal polynomials. b contains the coefficients 
% for the expansion  P(x) = b_1 p_1(x) + ... + b_m p_m(x).
% (Note that we used m=n+1)
m = numel(t);  % Maximal degree of orhtogonal polynomial
alpha(1) = sum(t)/m;
% Initialization of recursion
p1 = ones(size(t)); 
p2 = t-alpha(1);
% Store values of orthogonal polynomials 
P = [p1;p2];
% Compute first two expansion coefficients
b(1) = dot(p1,y)/dot(p1,p1);
b(2) = dot(p2,y)/dot(p2,p2);
% Main loop: 3-term recursion
for k=1:min(n-1,m-2)
  p0 = p1;  p1 = p2;
  alpha(k+1) =dot(p1,(t.*p1))/norm(p1)^2;
  beta(k) = (norm(p1)/norm(p0))^2;
  p2 = (t-alpha(k+1)).*p1-beta(k)*p0;
  % Store values
  P = [P; p2];
  % Compute next expansion coefficient
  b(k+2) = dot(p2,y)/norm(p2)^2;
end
