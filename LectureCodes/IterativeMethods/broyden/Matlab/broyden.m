function [x,prot] = broyden(F,x,J,tol,maxit,sol)
% Good Broyden rank-1-update quasi-Newton method
% Straightforward implementation for small problems
% F = Non-linear mapping in n dimensions
% x = initial guess
% J = initial guess for Jacobi matrix at x0
% tol = tolerance for termination

if (nargin < 6), sol = []; end
if (nargin < 5), maxit = 20; end

% Initial step 
k = 1; tol = tol^2; s = J\F(x); 
x = x - s; f = F(x); l = dot(s,s);

% Keeping a record of convergence history
prot = [k,sqrt(l),norm(f)];
if (~isempty(sol)), prot = [prot, norm(x-sol)]; end
while ((l > tol) && (k < 10))
  J = J - f*s'/l;
  s = J\f; x = x - s; 
  f = F(x); l = dot(s,s);
  fprintf('Iteration %d: |s| = %f, |F(x)| = %f \n',k,sqrt(l),norm(f));
  k = k+1;
  res = [k,sqrt(l),norm(f)];
  if (~isempty(sol)), res = [res, norm(x-sol)]; end
  prot = [prot; res];
end
