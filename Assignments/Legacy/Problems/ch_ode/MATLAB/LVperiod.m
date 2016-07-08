function [u0,v0]=LVperiod

%constants
alpha=2;
beta=1;
gamma=1;
delta=1;

% end time
T=5;

% initial condition
y_0=[3;2];

% right hand side
y_dot=@(t,y) [(alpha-beta*y(2))*y(1);(delta*y(1)-gamma)*y(2)];

% set the same tolerances as in exercise (6c)
opts=odeset('AbsTol',1e-8,'RelTol',1e-8);

% calculate the end point with the new initial value
% don't forget to set the right tolerances
[tout,yout]=ode45(y_dot,[0 T],y_0,opts);

% Jacobi matrix of the right hand side
DF=@(t,y) [(alpha-beta*y(2)), -beta*y(1);...
    delta*y(2), (delta*y(1)-gamma)];

% Newton's method

% termination criterion
while (norm(yout(end,:)'-y_0)>1e-5)
    
    % calculate Wronskian
    W=wronski(y_dot,DF,T,y_0);
    
    % carry out a step of Newton's method
    y_0=y_0-(W-eye(2))\(yout(end,:)'-y_0);
    
    % calculate the end point with the new intiial value
    % don't forget to set the right tolerances
    [tout,yout]=ode45(y_dot,[0 T],y_0,opts);
       
end

% save the result
u0=y_0(1);
v0=y_0(2);



end