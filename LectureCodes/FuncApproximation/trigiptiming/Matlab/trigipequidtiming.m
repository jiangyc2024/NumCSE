function trigipequidtiming
% Runtime comparison between efficient ($\to$ Code~\ref{trigipequid}) and direct computation
% ($\to$ Code~\ref{trigpolycoeff} of coefficients of trigonoetric interpolation polynomial in
% equidistant points.
Nruns = 3; times = [];
for n = 10:10:500
  disp(n)
  N = 2*n+1; t = 0:1/N:(1-1/N); y = exp(cos(2*pi*t));
  t1 = realmax; t2 = realmax; 
  for k=1:Nruns
    tic; [a,b] = trigpolycoeff(t,y); t1 = min(t1,toc);
    tic;  [a,b] = trigipequid(y); t2 = min(t2,toc);
  end
  times = [times; n , t1 , t2];
end
   
figure; loglog(times(:,1),times(:,2),'b+',...
               times(:,1),times(:,3),'r*');
xlabel('{\bf n}','fontsize',14);
ylabel('{\bf runtime[s]}','fontsize',14);
legend('trigpolycoeff','trigipequid','location','best');

print -depsc2 '../PICTURES/trigipequidtiming.eps';
