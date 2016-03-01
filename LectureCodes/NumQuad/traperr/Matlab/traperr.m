function traperr()

clear a;
global a;
l = 0; r = 0.5; % integration interval
N = 50;
a = 0.5; res05 = trapezoidal(@issin,l,r,1:N);
ex05 = trapezoidal(@issin,l,r,500); ex05 = ex05(1,2);
a = 0.9; res09 = trapezoidal(@issin,l,r,1:N);
ex09 = trapezoidal(@issin,l,r,500); ex09 = ex09(1,2);
a = 0.95; res95 = trapezoidal(@issin,l,r,1:N);
ex95 = trapezoidal(@issin,l,r,500); ex95 = ex95(1,2);
a = 0.99; res99 = trapezoidal(@issin,l,r,1:N);
ex99 = trapezoidal(@issin,l,r,500); ex99 = ex99(1,2);
figure('name','trapezoidal rule for non-periodic function');    
loglog(1./res05(:,1),abs(res05(:,2)-ex05),'r+-',...
     1./res09(:,1),abs(res09(:,2)-ex09),'b+-',...
     1./res95(:,1),abs(res95(:,2)-ex95),'c+-',...
     1./res99(:,1),abs(res99(:,2)-ex99),'m+-');
set(gca,'fontsize',12);
legend('a=0.5','a=0.9','a=0.95','a=0.99',3);
xlabel('{\bf no. of. quadrature nodes}','fontsize',14);
ylabel('{\bf |quadrature error|}','fontsize',14);
title('{\bf Trapezoidal rule quadrature for non-periodic function}','fontsize',12);

print -depsc2 '../PICTURES/traperr2.eps';

clear a;
global a;
l = 0; r = 1; % integration interval
N = 20;
a = 0.5; res05 = trapezoidal(@issin,l,r,1:N);
ex05 = trapezoidal(@issin,l,r,500); ex05 = ex05(1,2);
a = 0.9; res09 = trapezoidal(@issin,l,r,1:N);
ex09 = trapezoidal(@issin,l,r,500); ex09 = ex09(1,2);
a = 0.95; res95 = trapezoidal(@issin,l,r,1:N);
ex95 = trapezoidal(@issin,l,r,500); ex95 = ex95(1,2);
a = 0.99; res99 = trapezoidal(@issin,l,r,1:N);
ex99 = trapezoidal(@issin,l,r,500); ex99 = ex99(1,2);
figure('name','trapezoidal rule for periodic function');    
semilogy(1./res05(:,1),abs(res05(:,2)-ex05),'r+-',...
	 1./res09(:,1),abs(res09(:,2)-ex09),'b+-',...
	 1./res95(:,1),abs(res95(:,2)-ex95),'c+-',...
	 1./res99(:,1),abs(res99(:,2)-ex99),'m+-');
set(gca,'fontsize',12);
legend('a=0.5','a=0.9','a=0.95','a=0.99',3);
xlabel('{\bf no. of. quadrature nodes}','fontsize',14);
ylabel('{\bf |quadrature error|}','fontsize',14);
title('{\bf Trapezoidal rule quadrature for 1./sqrt(1-a*sin(2*pi*x+1))}','fontsize',12);

print -depsc2 '../PICTURES/traperr1.eps';

	 
