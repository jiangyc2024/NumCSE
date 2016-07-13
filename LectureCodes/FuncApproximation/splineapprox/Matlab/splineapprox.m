function splineapprox(f,df,a,b,N)
% Study of interpolation error norms for \emph{complete} cubic spline interpolation of \Blue{$f$} 
% on equidistant knot set. 
x = a:0.00025:b;  fv = feval(f,x); % fine mesh for norm evaluation
dfa = feval(df,a); dfb = feval(df,b); % Slopes at endpoints
err = [];
for j=2:N
  t = a:(b-a)/j:b;          % spline knots
  y = [dfa,feval(f,t),dfb];
  % compute complete spline imposing exact first derivative at the endpoints
  % Please refer to MATLAB documentation of \texttt{spline}
  v = spline(t,y,x);
  d = abs(fv-v);
  h = x(2:end)-x(1:end-1);
  % compute $L^2$ norm of the error using trapezoidal rule on mesh \texttt{x}
  l2 = sqrt(0.5*dot(h,(d(1:end-1).^2+d(2:end).^2)));
  % columns of err = meshwidth, $L^\infty$ error, $L^2$ error:
  err = [err; (b-a)/j,max(d),l2];
end

figure('Name','Spline interpolation');
plot(t,y(2:end-1),'m*',x,fv,'b-',x,v,'r-');
xlabel('{\bf t}','fontsize',14);  ylabel('{\bf s(t)}','fontsize',14);
legend('Data points','f','Cubic spline interpolant','location','best');

figure('Name','Spline approximation error');
loglog(err(:,1),err(:,2),'r.-',err(:,1),err(:,3),'b.-');
grid on;	xlabel('Meshwidth h');	  ylabel('||s-f||');
legend('sup-norm','L^2-norm', 'Location','NorthWest');
%compute algebraic orders of convergence using polynomial fit
p = polyfit(log(err(:,1)),log(err(:,2)),1); exp_rate_Linf=p(1)
p = polyfit(log(err(:,1)),log(err(:,3)),1); exp_rate_L2=p(1)

