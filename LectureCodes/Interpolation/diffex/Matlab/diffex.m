% 16.06.2009        diffex.m
% numeric differentiation through interpolation
% approximate the derivative of f in x0 with prescribed tolerance tol
% to be compared with evaldirnumdiff.m
% the elements of vector d are the approximations of f'(x_0)
%
% example:
% diffex(@tan,1.1,0.5, 10^(-10))



function d = diffex(f,x,h0,tol)
format long;
h = h0;
y = (feval(f,x+h0)-feval(f,x-h0))/(2*h0);
d = y;
for i=2:10
  h(i) = h(i-1)/2;
  y = [y,(feval(f,x+h(i))-feval(f,x-h(i)))/h(i-1)];
  for k=i-1:-1:1
    y(k) = y(k+1)-(y(k+1)-y(k))*h(i)/(h(i)-h(k));
  end
  d = [d;y(1)];
  if (abs(y(2)-y(1)) < tol*abs(y(1))), break; end
end
