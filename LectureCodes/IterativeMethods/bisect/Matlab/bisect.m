function x = bisect(@{\com{F}}@,a,b,tol)
% Searching zero of \Blue{$F$} in \Blue{$[a,b]$} by bisection
if (a>b), t=a; a=b; b=t; end;
fa = F(a); fb = F(b);
if (fa*fb>0), error('f(a), f(b) same sign'); end;
if (fa > 0), v=-1; else v = 1; end
x = 0.5*(b+a);
while((b-a > tol) && @\com{((a<x) \& (x<b))}@) @\Label[line]{bs:2}@
  if (v*F(x)>0), b=x; else a=x; end;
  x = 0.5*(a+b)
end  
