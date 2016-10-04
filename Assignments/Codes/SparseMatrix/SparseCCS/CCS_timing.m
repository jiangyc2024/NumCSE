nruns = 10; res = [];
ns = 2.^(5:12);
for n = ns
  A   = [n+1, (n:-1:1); ones(n,1), eye(n,n)];
  Asp = sparse(A);
  t  = realmax;
  t2 = realmax;
  for k=1:nruns
    tic; CCS(A);   t  = min(toc,t);
    tic; CCS(Asp); t2 = min(toc,t2);
  end
  res = [res; t t2];
end
figure('name','timings of CCS for full and sparse formats of matrix A');
c1 = sum(res(:,1))/sum(ns.^2);
c2 = sum(res(:,2))/sum(ns);
loglog(ns, res(:,1),'bo', ns, res(:,2), 'rx', ...
       ns, c1*ns.^2, 'g-', ns, c2*ns, 'k--');
xlabel('{\bf dimension n}','fontsize',14);
ylabel('{\bf time[s]}','fontsize',14);
title('{\bf timings of CCS for full and sparse formats of matrix A}','fontsize',14);
legend('full','sparse','O(n^2)','O(n)','location','best');
print -depsc2 '../PICTURES/CCS_timing.eps';
