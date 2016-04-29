function s=hermloceval(t,t1,t2,y1,y2,c1,c2)
% \Blue{$y_{1}$}, \Blue{$y_{2}$}: data values, \Blue{$c_{1}$}, \Blue{$c_{2}$}: slopes
h = t2-t1; t = (t-t1)/h;
a1 = y2-y1; a2 = a1-h*c1;
a3 = h*c2-a1-a2;
s = y1+(a1+(a2+a3*t).*(t-1)).*t;
