function hermiteapprox1(f,df,a,b,N)
% Investigation of interpolation error norms for cubic Hermite interpolation of \Blue{$f$} (handle \texttt{f})
% on \Blue{$[a,b]$} with slopes given by \Blue{$f'$} (function handle \texttt{df}).
% \texttt{N} gives the maximum number of mesh intervals
err = [];
for j=2:N
  xx=a; % \texttt{xx} is the fine mesh on which the error norms are computed
  val=f(a); % function values of mesh \texttt{xx}
  
  t = a:(b-a)/j:b; % mesh nodes
  y = f(t); c=df(t); % function values and slopes provide coefficients

  for k=1:j-1
    vx = linspace(t(k),t(k+1), 100);
    % See Code~\ref{hermloceval} for local evaluation of Hermite interpolant
    locval=hermloceval(vx,t(k),t(k+1),y(k),y(k+1),c(k),c(k+1));
    xx=[xx, vx(2:100)];
    val=[val,locval(2:100)];
  end
  d = abs(feval(f,xx)-val); % pointwise error on mesh \texttt{xx}
  h = (b-a)/j;
  % compute $L^2$ norm of the error using trapezoidal rule
  l2 = sqrt(h*(sum(d(2:end-1).^2)+(d(1)^2+d(end)^2)/2) );
  % columns of err = meshwidth, sup-norm error, $L^2$ error:
  err = [err; h, max(d),l2];
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

print -depsc2 '../PICTURES/hermiperrslopes.eps';