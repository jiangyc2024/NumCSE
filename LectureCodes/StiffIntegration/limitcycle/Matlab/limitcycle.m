fun = @(t,y) ([-y(2);y(1)] + lambda*(1-y(1)\symbol{94}2-y(2)\symbol{94}2)*y);
tspan = [0,2*pi];
y0 = [1,0];
opts = odeset('stats','on','reltol',1E-4,'abstol',1E-4);
[t45,y45] = ode45(fun,tspan,y0,opts);
