function LVtrajectories

%constants
alpha=2;
beta=1;
gamma=1;
delta=1;

%initial value 1
y_0=[4;2];

%right hand side
y_dot=@(t,y) [(alpha-beta*y(2))*y(1);(delta*y(1)-gamma)*y(2)];

%end time
T=5;

%integrate with ode45
[tout,yout]=ode45(y_dot,[0 T],y_0);

%plot the result
figure;
plot(yout(:,1),yout(:,2),'*-b'); %'b' stands for blue
hold on;
%initial value as circle
plot(yout(1,1),yout(1,2),'or'); %'or' stands for a red circle
%end point as square
plot(yout(end,1),yout(end,2),'sr'); %'s' stands for a square
%legend
legend('Initial value 1','Initial value 1','End point 1');

%initial value 2
y_0=[3;2];

%integrate with ode45
[tout,yout]=ode45(y_dot,[0 T],y_0);

%Plot the result
plot(yout(:,1),yout(:,2),'*-r');
%initial value as circle
plot(yout(1,1),yout(1,2),'ob');
%end point as square
plot(yout(end,1),yout(end,2),'sb');
%legend
legend('Initial value 1','Initial value 1','End point 1',...
    'Initial value 2','Initial value 2','End point 2');

end
