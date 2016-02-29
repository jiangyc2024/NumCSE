function res = trapezoidal(fnct,a,b,N)
% Numerical quadrature based on trapezoidal rule
% fnct handle to y = f(x)
% a,b bounds of integration interval
% N+1 = number of equidistant integration points (can be a vector)
res = [];
for n = N
  h = (b-a)/n;  x = (a:h:b); w = [0.5 ones(1,n-1) 0.5];
  res = [res; h, h*dot(w,feval(fnct,x))];
end
