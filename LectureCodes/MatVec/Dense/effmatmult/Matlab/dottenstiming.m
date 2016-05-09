function dottenstiming(N,nruns)
% This function compares the runtimes for the multiplication of a vector with a
% rank-1 matrix $\Va\Vb^{\top}$, $\Va,\Vb\in\bbR^n$ using different associative
% evaluations measurements consider minimal time for several (\texttt{nruns}) runs

if (nargin < 1), N = 2.^(2:13); end
if (nargin < 2), nruns = 10; end 

times = []; % matrix for storing recorded runtimes
for n=N
  % Initialize dense vectors $\Va$, $\Vb$, $\Vx$ (column vectors!)
  a = (1:n)'; b = (n:-1:1)'; x = rand(n,1);

  % Measuring times using MATLAB tic-toc commands
  tfool = realmax; for i=1:nruns,  tic; y = (a*b')*x; tfool = min(tfool,toc); end;  
  tsmart = realmax; for i=1:nruns,  tic; y = a*dot(b',x); tsmart = min(tsmart,toc); end;  
  times = [times;n, tfool, tsmart];
end

% log-scale plot for investigation of asymptotic complexity
figure('name','dottenstiming');
loglog(times(:,1),times(:,2),'m+',...
       times(:,1),times(:,3),'r*',...
       times(:,1),times(:,1)*times(1,3)/times(1,1),'k-',...
       times(:,1),(times(:,1).^2)*times(2,2)/(times(2,1)^2),'b--');
xlabel('{\bf problem size n}','fontsize',14);
ylabel('{\bf average runtime (s)}','fontsize',14);
title('tic-toc timing, mininum over 10 runs');
legend('slow evaluation','efficient evaluation',...
       'O(n)','O(n^2)','location','northwest');
print -depsc2 '../PICTURES/dottenstiming.eps';

  