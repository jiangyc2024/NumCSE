function [x,k] = fastbroyd(F,x,J,tol,maxit)
% Good Broyden rank-1-update quasi-Newton method
% Straightforward implementation for small problems
% F = Non-linear mapping in n dimensions
% x = initial guess
% J = initial guess for Jacobi matrix at x0
% tol = tolerance for termination

if (nargin < 5), maxit = 20; end

[L,U] = lu(J);

k = 1; tol = tol^2; s = U\(L\F(x));
x = x - s; f = F(x); sn = dot(s,s);

% Storing update vectors
dx = [s];
dxn = [sn];

while ((sn > tol) && (k < maxit))
  w = U\(L\f);
  for l=2:k
    w = w+dx(:,l)*(dx(:,l-1)'*w)/dxn(l-1);
  end
  z = s'*w; s = (1+z/(sn-z))*w; sn = s'*s; dx = [dx,s]; dxn = [dxn,sn];
  x = x - s; f = F(x); 
  k = k+1;
end
