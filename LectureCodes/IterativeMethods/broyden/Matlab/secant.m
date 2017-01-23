function x = secant(F,x,J,tol,maxit)
% Good Broyden rank-1-update quasi-Newton method
% F = Non-linear mapping in n dimensions
% x = initial guess
% J = initial guess for Jacobi matrix at x0
% tol = tolerance for termination

if (nargin < 5), maxit = 20; end

k = 1;
s = J\F(x); x1 = x - s; 
fprintf('SV-Iteration 0: |s| = %f, |F(x)| = %f, x= %f, J = %f\n',...
	abs(s),abs(F(x1)),x1,J);
while ((norm(s) > tol) && (k < 10))
  J = (F(x1)-F(x))/(x1-x);
  s = F(x1)/J;
  x = x1; x1 = x-s;
  fprintf('SV-Iteration %d: |s| = %f, |F(x)| = %f, x= %f, J = %f\n',...
	  k,abs(s),abs(F(x1)),x,J);
  k = k+1;
end
