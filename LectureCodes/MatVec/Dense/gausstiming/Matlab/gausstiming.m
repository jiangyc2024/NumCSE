% MATLAB script for timing numerical solution of linear systems
nruns = 3; times = [];
for n = 2.^(3:12)
  fprintf('Matrix size n = %d\n',n);
  % Initialized random matrix and right hand side
  A = rand(n,n) + n*eye(n); b = rand(n,1);
  t1 = realmax; t2 = realmax; 
  for j=1:nruns
    tic; x1 = gausselimsolve(A,b); t1 = min(t1,toc);
    tic; x2 = A\b; t2 = min(t2,toc);
      norm(x1-x2),
  end
  times = [times; n t1 t2];  
end
    
figure('name','gausstiming');
loglog(times(:,1),times(:,2),'r-+',times(:,1),times(:,3),'m-*',...
       times(:,1),times(:,1).^3*(1E-5/(times(1,1)^3)),'k-');
xlabel('matrix size n','fontsize',14);
ylabel('execution time [s]','fontsize',14);
legend('gausselimsolve','backslash solve','O(n^3)','location','northwest');

print -depsc2 '../PICTURES/gausstiming.eps';
