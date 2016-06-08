%A = Wilkinsons counter-example for stable pivot strategy
%c = 2-condition of matrix A
%L, U, P = lower, upper and permutation matrix of lu-decomposition

%There seems to be a problem with Matlab ver. 7.0.1 in 
%combination with Intel Pentium IV processors. Somehow there's
%an overflow or something, big errors are made!

res = [];
for n=10:10:1000
  A = [tril(-ones(n,n-1))+2*[eye(n-1); zeros(1,n-1)],ones(n,1)];
  c = cond(A,2);
  [L,U,P] = lu(A);
  E = L*U-A;
  e = norm(E,2)/norm(A,2);
  res = [res; n,c,e];
end

figure('name','Wilkinson example: condition number');
plot(res(:,1),res(:,2),'r-+')
xlabel('n','Fontsize',14);
ylabel('cond_2(A)','Fontsize',16);

print -depsc2 '../PICTURES/wilkcond.eps';

figure('name','Wilkinson example: error in LU-decomposition');
plot(res(:,1),res(:,3),'m-+');
xlabel('{\bf n}','Fontsize',14);
ylabel('E','Fontsize',16);
legend(sprintf('On Intel Core i7, Matlab ver %s',version),'location','northwest');

print -depsc2 '../PICTURES/wilkerr.eps';
