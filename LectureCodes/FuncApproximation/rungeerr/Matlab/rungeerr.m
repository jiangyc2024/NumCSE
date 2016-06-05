% interpolation error plot for Runge's example
% (``Quick and dirty'' MATLAB implementation, see \ref{subsec:interp-algorithms})
f = @(x) (1./(1+x.^2));
x = -5:0.01:5;  % sampling point for approximate maximum norm
fv = f(x);
err = [];
for d=1:20
  t = -5+(0:d)*10/d;
  p = polyfit(t,f(t),d);
  y = polyval(p,x);
  err = [err; d, max(abs(y-fv))];
end

figure; semilogy(err(:,1),err(:,2),'r--+');
xlabel('{\bf degree d}','fontsize',14);
ylabel('{\bf interpolation error (maximum norm)}','fontsize',14);
print -depsc2 '../PICTURES/rungeerrmax.eps';

