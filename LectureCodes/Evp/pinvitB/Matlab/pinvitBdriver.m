function pinvitBdriver

k = 1; tol = 0; maxit = 50;
for n=[50,100,200]
  A = spdiags(repmat([1/n,-1,2*(1+1/n),-1,1/n],n,1),[-n/2,-1,0,1,n/2],n,n);
  evalA = @(x) A*x;
  lmax = min(eig(full(A)))
  
%   res = []; z = rand(n,1); 
%   for l=1:150, z = A\z; z = z/norm(z); res = [res;dot(A*z,z) - lmax]; end
%   res,z

% inverse iteration
  invB = @(x) A\x;
  [lmin,z,rn{k}] = pinvitB(evalA,n,invB,tol,maxit);
  nit{k} = size(rn{k},1); err{k} = abs(rn{k} - lmax)/lmax;
  
% tridiagonal preconditioning
  B = spdiags(spdiags(A,[-1,0,1]),[-1,0,1],n,n);
  invB  = @(x) B\x;
  [lmin,z,rnpc{k}] = pinvitB(evalA,n,invB,tol,maxit);
  nitpc{k} = size(rnpc{k},1); errpc{k} = abs(rnpc{k} - lmax)/lmax;
  k = k+1;
end 

figure('name','error norms');
semilogy(0:nit{1}-1,err{1},'+b');   hold on;
semilogy(0:nit{2}-1,err{2},'+m');
semilogy(0:nit{3}-1,err{3},'+r');
semilogy(0:nitpc{1}-1,errpc{1},'*b');
semilogy(0:nitpc{2}-1,errpc{2},'*m');
semilogy(0:nitpc{3}-1,errpc{3},'*r');
legend('INVIT, n = 50','INVIT, n = 100','INVIT, n = 200',...
       'PINVIT, n = 50','PINVIT, n = 100','PINVIT, n = 200','location','northeast');
xlabel('{\bf # iterationstep}','fontsize',14);
ylabel('{\bf error in approximation for \lambda_{max}}','fontsize',14);

print -depsc2 '../PICTURES/pinvitB.eps';

tol = 1E-4;
itnum = [];
for n = 2.^(4:15)
  A = spdiags(repmat([1/n,-1,2*(1+1/n),-1,1/n],n,1),[-n/2,-1,0,1,n/2],n,n);
  evalA = @(x) A*x;
  
% inverse iteration
  invB = @(x) A\x;
  [lmin,z,rn] = pinvitB(evalA,n,invB,tol,maxit);
  
% tridiagonal preconditioning
  B = spdiags(spdiags(A,[-1,0,1]),[-1,0,1],n,n);
  invB  = @(x) B\x;
  [lmin,z,rnpc] = pinvitB(evalA,n,invB,tol,maxit);
  
  itnum = [itnum; n, size(rn,1), size(rnpc,1)];
end

figure('name','iteration numbers');
semilogx(itnum(:,1),itnum(:,2),'b+',itnum(:,1),itnum(:,3),'r+');
xlabel('{\bf n}','fontsize',14);
ylabel('{\bf #iteration steps}','fontsize',14);
title('(P)INVIT iterations: tolerance = 0.0001');
legend('INVIT','PINVIT','location','northwest');

print -depsc2 '../PICTURES/pinvitBitnum.eps';
