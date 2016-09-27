% Timing for row/column operations on matrices
% We conduct K runs in order to reduce the risk of skewed maesurements
% due to OS activity during MATLAB run.
K = 3; res = [];
for n=2.^(5:13)
  A = randn(n,n);
  
  t1 = realmax; 
  for k=1:K,  tic; 
    for j = 1:n-1,  A(:,j+1) = A(:,j+1) - A(:,j); end; 
    t1 = min(toc,t1); 
  end
  t2 = realmax; 
  for k=1:K,  tic; 
    for i = 1:n-1,  A(i+1,:) = A(i+1,:) - A(i,:); end;
    t2 = min(toc,t2); 
  end
  res = [res; n, t1 , t2];
end
save('../Eigen/timing_matlab.dat','res','-ascii')

% Plot runtimes versus matrix sizes
% figure; plot(res(:,1),res(:,2),'r+', res(:,1),res(:,3),'m*');
% xlabel('{\bf n}','fontsize',14); 
% ylabel('{\bf runtime [s]}','fontsize',14); 
% legend('A(:,j+1) = A(:,j+1) - A(:,j)','A(i+1,:) = A(i+1,:) - A(i,:)',...
%        'location','northwest');
% print -depsc2 '../PICTURES/accessrtlin.eps';
% 
% figure; loglog(res(:,1),res(:,2),'r+', res(:,1),res(:,3),'m*');
% xlabel('{\bf n}','fontsize',14); 
% ylabel('{\bf runtime [s]}','fontsize',14); 
% legend('A(:,j+1) = A(:,j+1) - A(:,j)','A(i+1,:) = A(i+1,:) - A(i,:)',...
%        'location','northwest');
% print -depsc2 '../PICTURES/accessrtlog.eps';
