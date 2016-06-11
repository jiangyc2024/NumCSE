function eigtiming

A = rand(500,500); B = A'*A; 
C = gallery('tridiag',500,1,3,1);
times = [];
for n=5:5:500
  An = A(1:n,1:n); Bn = B(1:n,1:n); Cn = C(1:n,1:n);
  t1 = 1000; for k=1:3, tic; d = eig(An); t1 = min(t1,toc); end
  t2 = 1000; for k=1:3, tic; [V,D] = eig(An); t2 = min(t2,toc); end
  t3 = 1000; for k=1:3, tic; d = eig(Bn); t3 = min(t3,toc); end
  t4 = 1000; for k=1:3, tic; [V,D] = eig(Bn); t4 = min(t4,toc); end
  t5 = 1000; for k=1:3, tic; d = eig(Cn); t5 = min(t5,toc); end
  times = [times; n t1 t2 t3 t4 t5];
end

figure;
loglog(times(:,1),times(:,2),'r+', times(:,1),times(:,3),'m*',...
       times(:,1),times(:,4),'cp', times(:,1),times(:,5),'b^',... 
       times(:,1),times(:,6),'k.');
xlabel('{\bf matrix size n}','fontsize',14);
ylabel('{\bf time [s]}','fontsize',14);
title('eig runtimes');
legend('d = eig(A)','[V,D] = eig(A)','d = eig(B)','[V,D] = eig(B)','d = eig(C)',...
       'location','northwest');

print -depsc2 '../PICTURES/eigtimingall.eps'

figure;
loglog(times(:,1),times(:,2),'r+', times(:,1),times(:,3),'m*',...
	 times(:,1),(times(:,1).^3)/(times(1,1)^3)*times(1,2),'k-');
xlabel('{\bf matrix size n}','fontsize',14);
ylabel('{\bf time [s]}','fontsize',14);
title('nxn random matrix');
legend('d = eig(A)','[V,D] = eig(A)','O(n^3)','location','northwest');

print -depsc2 '../PICTURES/eigtimingA.eps'

figure;
loglog(times(:,1),times(:,4),'r+', times(:,1),times(:,5),'m*',... 
	 times(:,1),(times(:,1).^3)/(times(1,1)^3)*times(1,2),'k-');
xlabel('{\bf matrix size n}','fontsize',14);
ylabel('{\bf time [s]}','fontsize',14);
title('nxn random Hermitian matrix');
legend('d = eig(A)','[V,D] = eig(A)','O(n^3)','location','northwest');	 

print -depsc2 '../PICTURES/eigtimingB.eps'
	 
figure;
loglog(times(:,1),times(:,6),'r*',...
	 times(:,1),(times(:,1).^2)/(times(1,1)^2)*times(1,2),'k-');
xlabel('{\bf matrix size n}','fontsize',14);
ylabel('{\bf time [s]}','fontsize',14);
title('nxn tridiagonel Hermitian matrix');
legend('d = eig(A)','O(n^2)','location','northwest');

print -depsc2 '../PICTURES/eigtimingC.eps'


	 
