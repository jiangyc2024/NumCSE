function hermintp1(f,t)
% compute and plot the cubic Hermite interpolant of the function \texttt{f} in the nodes \texttt{t}
% using weighted averages according to \eqref{pwintp:AverageSlopes} as local slopes
n = length(t); h = diff(t); % computes lengths of intervals between nodes: \Blue{$h_{i} := t_{i+1}-t_{i}$}  
y = f(t); % \texttt{f} must support collective evaluation for row vector argument
delta = diff(y)./h; % slopes of piecewise linear interpolant
c = [delta(1),...                          
     ((h(2:end).*delta(1:end-1)+h(1:end-1).*delta(2:end))...
        ./(t(3:end) - t(1:end-2)) ),...  
     delta(end)]; % slopes from weighted average, see \eqref{pwintp:AverageSlopes}

figure('Name','Hermite Interpolation');
plot(t,y,'ko'); hold on;  %  plot data points
fplot(f,[t(1), t(n)]);
for j=1:n-1  % compute and plot the Hermite interpolant with slopes \texttt{c}
  vx = linspace(t(j),t(j+1), 100);
  plot(vx,hermloceval(vx,t(j),t(j+1),y(j),y(j+1),c(j),c(j+1)),'r-', 'LineWidth',2);
end
for j=2:n-1  % plot segments indicating the slopes \Blue{$c_i$}
  plot([t(j)-0.3*h(j-1),t(j)+0.3*h(j)],...
       [y(j)-0.3*h(j-1)*c(j),y(j)+0.3*h(j)*c(j)],'k-','LineWidth',2);
end
plot([t(1),t(1)+0.3*h(1)],[y(1),y(1)+0.3*h(1)*c(1)],'k-', 'LineWidth',2);
plot([t(end)-0.3*h(end),t(end)],[y(end)-0.3*h(end)*c(end),y(end)],'k-', 'LineWidth',2);
xlabel('t'); ylabel('s(t)');
legend('Data points','f(x)','Piecew. cubic interpolation polynomial ');
hold off;
