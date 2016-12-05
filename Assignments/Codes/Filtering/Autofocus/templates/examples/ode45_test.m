% T = 10.;
% a1 = 3;
% a2 = 2;
% b1 = 0.1;
% b2 = 0.1;
% f = @(t, y) [(a1 - b1*y(2))*y(1); (b2*y(1) - a2)*y(2)];
% y0 = [100., 5.];
% tspan = [0., T];
% 
% tic
% [t45, y45] = ode45(f, tspan, y0);
% toc
% 
% plot(t45, y45);
% disp length(t45)

T = 10000.;
f = @(t, y) 1 / y;
y0 = 0.2;
tspan = [0., T];

options = odeset('AbsTol',1e-8,'RelTol',1e-6,'NormControl','on');
tic
[t45, y45] = ode45(f, tspan, y0, options);
toc

t45(end)
y45(end)
plot(t45, y45);
length(t45)