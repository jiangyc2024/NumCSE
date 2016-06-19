figure; spy(spdiags(repmat([-1 2 5],16,1),[-8,0,8],16,16)); % \label{spf:1}
title('Pattern for matrix {\bf A} for n = 16','fontsize',14);
print -depsc2 '../PICTURES/spdiagsmatspy.eps';

t = [];
for i=1:20
  n = 2^i; m = n/2;
  A = spdiags(repmat([-1 2 5],n,1),[-n/2,0,n/2],n,n); % \label{spf:2}
  
  t1 = inf; for j=1:5, tic; v = A(m,:)+j; t1 = min(t1,toc); end
  t2 = inf; for j=1:5, tic; v = A(:,m)+j; t2 = min(t2,toc); end
  t = [t; size(A,1), nnz(A), t1, t2 ];
end

figure;
loglog(t(:,1),t(:,3),'r+-', t(:,1),t(:,4),'b*-',...
       t(:,1),t(1,3)*t(:,1)/t(1,1),'k-');
xlabel('{\bf size n of sparse quadratic matrix}','fontsize',14);
ylabel('{\bf access time [s]}','fontsize',14);
legend('row access','column access','O(n)','location','northwest');

print -depsc2 '../PICTURES/sparseaccess.eps';

