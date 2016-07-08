function W=wronski(f,DF,T,y0)

% it will be useful to save the length of y0,
% as we will need it multiple times
n=length(y0);

%initial condition of the symstem of differential equation+variational equation
% the intiial value of the variational equation is the unity matrix, which
% must be transformed to a vector using the command 'reshape'
y0=[y0;reshape(eye(2),n^2,1)];

% right hand side of the system differential eq.+variational eq.
% the first 'n' components are the right hand side fo the diff'eq.
% About the right hand side of the variational eq.: firt transform the vector
% to a matrix, then multiply by 'DF(t,y(...))' and then
% transform the result back to a vector
fun=@(t,y) [f(t,y(1:n));...
    reshape(DF(t,y(1:n))*reshape(y(n+1:end),n,n),n^2,1)];

% set low tolerances
opts=odeset('AbsTol',1e-8,'RelTol',1e-8);

%integrate with ode45
[tout,yout]=ode45(fun,[0 T],y0,opts);

%read out Wronski-Matrix from yout,
%we have to apply 'reshape' once more
W=reshape(yout(end,n+1:end)',n,n);


end