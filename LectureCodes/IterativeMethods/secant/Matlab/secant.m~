function x = secant(x0,x1,F,rtol,atol,MAXIT)
fo = F(x0);
for i=1:MAXIT
  fn = F(x1);
  s = fn*(x1-x0)/(fn-fo); % correction
  x0 = x1; x1 = x1-s;
  if ((abs(s) < max(atol,rtol*min(abs([x0;x1])))))
    x = x1; return; end
  fo = fn; 
end
