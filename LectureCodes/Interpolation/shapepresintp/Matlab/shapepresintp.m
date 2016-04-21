% 22.06.2009            shapepresintp.m
% Shape preserving interpolation through nodes (t,y)
% Build a quadratic spline interpolation without overshooting
% In 4 steps, see the comments in each step
%
% example:
% shapepresintp([1, 1.5,3:10], [ 1 2.5  4 3.8  3 2 1  4 2 1])
% shapepresintp([0:0.05:1], cosf([0:0.05:1]))
% shapepresintp([0:0.05:1], [zeros(1,5),1,zeros(1,15) ] )
% shapepresintp([0:12], [0 0.3 0.5 0.2 0.6 1.2 1.3 1 1 1 1 0 -1])


function shapepresintp(t,y)

n=length(t)-1;
% plot the data and prepare the figure
figure;
hplot(1)=plot(t,y,'k*-');
hold on;
newaxis=axis;
newaxis([3,4])=[newaxis(3)-0.5, newaxis(4)+0.5];
axis(newaxis);          % enlarge the vertical size of the plot
title('Data points - Press enter to continue')
plot(t,ones(1,n+1)* (newaxis(3)+0.25),'k.');
set(gca, 'XTick', t)
leg={'Data','Slopes','Middle points','Linear spline' ,'Sh. pres. spline'};
legend(hplot(1),leg{1});
pause;

% ===== Step 1: choice of slopes =====
% shape-faithful slopes (c) in the nodes using harmonic mean of data slopes
% the final interpolant will be tangents to these segments

disp('STEP 1')
title('Shape-faithful slopes - Press enter to continue')
h=diff(t);
delta = diff(y) ./ h;             % slopes of data
c=zeros(size(t));
for j=1:n-1
  if (delta(j)*delta(j+1) >0)
    c(j+1) = 2/(1/delta(j) + 1/delta(j+1));
  end
end
c(1)=2*delta(1)-c(2);  c(n+1)=2*delta(n)-c(n);

% plot segments indicating the slopes c(i):
% use (vector) plot handle 'hplot' to reduce the linewidth in step 2
hplots=zeros(1,n+1);
for j=2:n
  hplotsl(j)=plot([t(j)-0.3*h(j-1),t(j)+0.3*h(j)], [y(j)-0.3*h(j-1)*c(j),y(j)+0.3*h(j)*c(j)],'-','linewidth',2);
end
hplotsl(1)=plot([t(1),t(1)+0.3*h(1)], [y(1),y(1)+0.3*h(1)*c(1)], '-','linewidth',2);
hplotsl(n+1)= plot([t(end)-0.3*h(end),t(end)], [y(end)-0.3*h(end)*c(end),y(end)], '-','linewidth',2);
legend([hplot(1), hplotsl(1)], leg{1:2});
pause;

% ===== Step 2: choice of middle points =====
% fix points p(j) in [t(j), t(j+1)], depending on the slopes c(j), c(j+1)

disp('STEP 2')
title('Middle points - Press enter to continue')
set(hplotsl,'linewidth',1)

p = (t(1)-1)*ones(1,length(t)-1);
for j=1:n
  if (c(j) ~= c(j+1))
    p(j)=(y(j+1)-y(j)+ t(j)*c(j)-t(j+1)*c(j+1)) / (c(j)-c(j+1));
  end
  % check and repair if p(j) is outside its interval:
  if ((p(j)<t(j))||(p(j)>t(j+1)));    p(j) = 0.5*(t(j)+t(j+1));  end;
end

hplot(2)=plot(p,ones(1,n)*(newaxis(3)+0.25),'go');
legend([hplot(1), hplotsl(1), hplot(2)], leg{1:3});
pause;

% ===== Step 3: auxiliary linear spline =====
% build the linear spline with nodes in:
% -t(j)
% -the middle points between t(j) and p(j) 
% -the middle points between p(j) and t(j+1) 
% -t(j+1) 
% and with slopes c(j) in t(j), for every j 

disp('STEP 3')
title('Auxiliary linear spline - Press enter to continue')

for j=1:n
  hplot(3)=plot([t(j) 0.5*(p(j)+t(j)) 0.5*(p(j)+t(j+1)) t(j+1)],    [y(j) y(j)+0.5*(p(j)-t(j))*c(j) y(j+1)+0.5*(p(j)-t(j+1))*c(j+1) y(j+1)],  'm-^');
  plot ([t(j) 0.5*(p(j)+t(j)) 0.5*(p(j)+t(j+1)) t(j+1)],  ones(1,4)*(newaxis(3)+0.25)  ,  'm^');
end
legend([hplot(1), hplotsl(1), hplot(2), hplot(3)], leg{1:4});
pause;


% ===== Step 4: quadratic spline =====
% final quadratic shape preserving spline
% quadratic polynomial in the intervals  [t(j), p(j)] and  [p(j), t(j)]
% tangent in t(j) and p(j) to the linear spline of step 3

disp('STEP 4')
title('Quadratic spline')

% for every interval 2 quadratic interpolations
% a, b, ya, yb = extremes and values in the subinterval
% w = value in middle point that gives the right slope
for j=1:n
    a=t(j);
    b=p(j);
    ya = y(j); 
    w = y(j)+0.5*(p(j)-t(j))*c(j);
    yb = ((t(j+1)-p(j))*(y(j)+0.5*(p(j)-t(j))*c(j))+...
      (p(j)-t(j))*(y(j+1)+0.5*(p(j)-t(j+1))*c(j+1)))/(t(j+1)-t(j));
    x=linspace(a,b,100);
    pb = (ya*(b-x).^2 + 2*w*(x-a).*(b-x)+yb*(x-a).^2)/((b-a)^2);
    hplot(4)=plot(x,pb,'r-', 'linewidth',2);  
  
    a = b;
    b = t(j+1);
    ya = yb;
    yb = y(j+1);
    w = y(j+1)+0.5*(p(j)-t(j+1))*c(j+1);
    x = (a:(b-a)/100:b);
    pb = (ya*(b-x).^2 + 2*w*(x-a).*(b-x)+yb*(x-a).^2)/((b-a)^2);
    plot(x,pb,'r-', 'linewidth',2);
    
    plot(p(j),ya,'go');
end

% replot initial nodes over the other plots:
plot(t,y,'k*');                                 
% plot(p,yb,'go')
legend([hplot(1), hplotsl(1), hplot(2:4)], leg);
title('Shape preserving interpolation')

