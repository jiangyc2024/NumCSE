function [x,k] = newtontol(x,F,DF,tol,maxit)
% MATLAB template for Newton method
s = feval(DF,x)\feval(F,x); x = x-s; k = 1;
while ((norm(s) > tol) && (k < maxit))
  s = DF(x)\F(x); x = x-s; k = k+1;
end
