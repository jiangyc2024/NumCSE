function Ival = adaptquad_nc(f,x,reltol,abstol)
global res;

m = length(x)-1; if (m > res.mmax), error('mmax reached'); end

h = diff(x); mp = 0.5*(x(1:end-1)+x(2:end));
fx = f(x); fm = f(mp);

trp_loc  = 0.5*h.*(0.5*fx(1:end-1)+fm+0.5*fx(2:end));
simp_loc = h.*(fx(1:end-1)+4*fm+fx(2:end))/6;
Ival = sum(simp_loc);
est_loc  = abs(simp_loc -trp_loc);
err_tot  = sum(est_loc);

pause;
res.level = res.level + 1;
res.X{res.level} = x;
res.vals(res.level) = Ival;
res.ests(res.level) = err_tot;
fprintf('** level %d: ',res.level);
fprintf('%d points, Ival = %f, err = %f\n',m,Ival,err_tot);
N = 1000; xp = x(1):(x(end)-x(1))/N:x(end);
[AX,H1,H2] = plotyy(xp,f(xp),mp,est_loc/max(est_loc),'plot','bar');
set(H1,'color',[0 1 0],'linewidth',2);
set(H2,'facecolor',[1 0 0]);
set(get(AX(1),'Ylabel'),'String','{\bf f}');
set(get(AX(1),'Xlabel'),'String','{\bf x}');
set(get(AX(2),'Ylabel'),'String','{\bf normalized estimated local error}');

if ((err_tot > reltol*abs(Ival)) && (err_tot > abstol))
  refcells = find(est_loc > 0.9*sum(est_loc)/length(h));
  Ival = adaptquad_nc(f,sort([x,mp(refcells)]),reltol,abstol);
end


