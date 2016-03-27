fun = @(t,y) y.^2;
[t1,y1] = ode45(fun,[0 2],1);
[t2,y2] = ode45(fun,[0 2],0.5);
[t3,y3] = ode45(fun,[0 2],2);
