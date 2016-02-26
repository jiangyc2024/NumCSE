function x = sqrtit(a)
x_old = -1; x = a;
while (x_old ~= x)
  x_old = x;
  x = 0.5*(x+a/x);
end
