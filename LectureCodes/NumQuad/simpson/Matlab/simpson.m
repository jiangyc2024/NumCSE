function res = simpson(fnct,a,b,N)
% Numerical quadrature based on Simpson rule
% fnct handle to y = f(x)
% a,b bounds of integration interval
% N+1 = number of equidistant integration points (can be a vector)

res = [];
for n = N
  h = (b-a)/n;
  x = (a:h/2:b);
  fv = feval(fnct,x);
  val = sum(h*(fv(1:2:end-2)+4*fv(2:2:end-1)+fv(3:2:end)))/6;
  res = [res; h, val];
end
