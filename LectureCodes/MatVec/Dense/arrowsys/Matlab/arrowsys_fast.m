function x = arrowsys_fast(d,c,b,alpha,y)
z = c./d;          % \Blue{$\Vz = \VD^{-1}\Vc$}  
w = y(1:end-1)./d; % \Blue{$\Vw = \VD^{-1}\Vy_1$}
xi = (y(end)-dot(b,w))/(alpha - dot(b,z));
x = [ w-xi*c./d;xi];
