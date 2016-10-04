nruns = 3; res = [];
ns = 2.^(2:12);
for n = ns
  d1 = (1:n)'; d2 = -d1; c = ones(n,1); b = [d1;d1];
  t  = realmax;
  for k=1:nruns
    tic; x  = solveA(d1,c,d2,b); t  = min(toc,t);
  end
  res = [res; t];
end
figure('name','timings for solveA');
c1 = sum(res)/sum(ns);
loglog(ns, res,'bo', ns, c1*ns, 'g-');
xlabel('{\bf dimension n}','fontsize',14);
ylabel('{\bf time[s]}','fontsize',14);
title('{\bf timings for solveA}','fontsize',14);
legend('solveA','O(n)','location','best');
print -depsc2 '../PICTURES/solveA_timing.eps';
