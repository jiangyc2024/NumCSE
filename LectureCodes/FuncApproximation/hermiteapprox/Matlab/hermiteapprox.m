% 16.11.2009              hermiteapprox.m
% Plot convergence of approximation error of cubic hermite interpolation
% with respect to the meshwidth
% print the algebraic order of convergence in sup and $L^2$ norms
% Slopes: weighted average of local slopes
%
% inputs:   f       function to be interpolated
%           a, b    left and right extremes of the interval
%           N       maximum number of subdomains

function hermiteapprox(f,a,b,N)
% Investigation of interpolation error norms for cubic Hermite interpolation of \Blue{$f$} (handle \texttt{f})
% on \Blue{$[a,b]$} with linearly averaged slopes according to \eqref{pwintp:AverageSlopes}.
% \texttt{N} gives the maximum number of mesh intervals
err = [];
for j=2:N
  xx=a;     % \texttt{xx} is the fine mesh on which the error norms are computed
  val=f(a); % function values on \texttt{xx}
  
  t = a:(b-a)/j:b;     % mesh nodes
  y = f(t); c=slopes1(t,y); % coefficients for Hermit polynomial representation

  for k=1:j-1
    vx = linspace(t(k),t(k+1), 100);
    locval=hermloceval(vx,t(k),t(k+1),y(k),y(k+1),c(k),c(k+1));
    xx=[xx, vx(2:100)];
    val=[val,locval(2:100)];
  end
  d = abs(feval(f,xx)- val);
  h = (b-a)/j;
  % compute $L^2$ norm of the error using trapezoidal rule
  l2 = sqrt(h*(sum(d(2:end-1).^2)+(d(1)^2+d(end)^2)/2) ); 
  % columns of err = meshwidth, sup-norm error, $L^2$ error:
  err = [err; h,max(d),l2];
end

figure('Name','Hermite approximation');
loglog(err(:,1),err(:,2),'r.',err(:,1),err(:,3),'b.');
grid on;
xlabel('{\bf meshwidth h}','fontsize',14);
ylabel('{\bf norm of interpolation error}','fontsize',14);
legend('sup-norm','L^2-norm','location','northwest');

% compute estimates for algebraic orders of convergence
% using linear regression on half of the data points
pI = polyfit(log(err(ceil(N/2):N-2,1)),log(err(ceil(N/2):N-2,2)),1);  exp_rate_Linf=pI(1)
pL2 = polyfit(log(err(ceil(N/2):N-2,1)),log(err(ceil(N/2):N-2,3)),1); exp_rate_L2=pL2(1)
hold on;
plot([err(1,1),err(N-1,1)], [err(1,1),err(N-1,1)].^pI(1)*exp(pI(2)),'m')
plot([err(1,1),err(N-1,1)], [err(1,1),err(N-1,1)].^pL2(1)*exp(pL2(2)),'k')

print -depsc2 '../PICTURES/hermiperravgsl.eps';

%---------------------------------------------
function    c=slopes1(t,y)
h = diff(t);                                % increments in t
delta = diff(y)./h;                         % slopes of piecewise linear interpolant
c = [delta(1),...                          
     ((h(2:end).*delta(1:end-1)+h(1:end-1).*delta(2:end))...
        ./(t(3:end) - t(1:end-2)) ),...
     delta(end)];



