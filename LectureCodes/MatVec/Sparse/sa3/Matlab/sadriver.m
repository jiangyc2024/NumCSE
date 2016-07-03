% compares running time for solving
% a problem with arrow matrices A:
% x = A\b

K = 3;
t = [];
for i=3:12
  n = 2^i; alpha = 2;
  b = ones(n,1); c = (1:n)';
  d = -ones(n,1); y = (-1).^(1:(n+1))';
  t1 = 1000; 
  for k=1:K
    tic;
    x1 = solvearrow1(alpha,b,c,d,y);
% naive solving 
    t1 = min(t1,toc);
  end
  t2 = 1000;
  for k=1:K
    tic;
    x2 = solvearrow2(alpha,b,c,d,y);
% structure sensible solving
    t2 = min(t2,toc);
  end
  t3 = 1000;
  for k=1:K
    tic;
    x3 = solvearrow3(alpha,b,c,d,y);
% solving with sparse matrices
    t3 = min(t3,toc);
  end
  t = [t; n t1 t2 t3];
  max(abs(x1-x3))
end

figure;
loglog(t(:,1),t(:,2),'b-*',t(:,1),t(:,3),'r-+');
xlabel('{\bf n}','Fontsize',14);
ylabel('{\bf time[seconds]}');
legend('Implementation I','Implementation II','location','northwest');

print -depsc2 '../PICTURES/satimes.eps';

figure;
loglog(t(:,1),t(:,2),'b-*',t(:,1),t(:,3),'r-+',t(:,1),t(:,4),'m-+');
xlabel('n','Fontsize',14);
ylabel('time[seconds]');
legend('Implementation I','Implementation II','Sparse Solver','location','northwest');

print -depsc2 '../PICTURES/sasp.eps';
