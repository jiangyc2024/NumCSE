% Measuring time

ksteps = 4;
Nmax = 1500;
prot = [];
for n=10:10:Nmax
  fprintf('<<<<< n=%d >>>>>>>>\n',n);
  s = ones(n,1); % Exact solution
  b = (1:n)';
  y = b-s;
  a = 1/sqrt(sum(y))*y;
  A = (eye(n) + a*a');
  F = @(x) (diag(x)*A*x - b);
  DF = @(x) (diag(x)*A + diag(A*x));
% Initial guess
  a = 2; b = 4; h = (b-a)/n;
  x0 = (a:h:b-h)';
  tol = 2*n*(1.0e-5);

  disp('Measuring time for Broydens method');
  bt = realmax;
  for k=1:ksteps;
    tic;
    [x,bk] = fastbroyd(F,x0,DF(x0),tol,20);
    t = toc;
    fprintf('Broyden(n=%d): %d steps, time = %e, final residual F(s) = %e\n',...
	    n,bk,t,norm(F(x)));
    bt = min(t,bt);
  end
  
  disp('Measuring time for Newtons method');
  nt = realmax;
  for k=1:ksteps;
    tic;
    [x,nk] = newtontol(x0,F,DF,tol,20);
    t = toc;
    fprintf('Newton(n=%d): %d steps, time = %e, final residual F(s) = %e\n',...
	    n,nk,t,norm(F(x)));
    nt = min(t,nt);
  end
  prot = [prot;n,bk,bt,nk,nt];
end

figure('name','Timing: #(steps)');
plot(prot(:,1),prot(:,2),'rs',prot(:,1),prot(:,4),'bs');
axis([0 Nmax+10 0 20]); 
set(gca,'fontsize',14);
grid on;
xlabel('n');
ylabel('Anzahl Schritte');
legend('Broyden-Verfahren','Newton-Verfahren');

print -dpsc2 'bigbroydensteps.eps';

figure('name','Timing: times');
plot(prot(:,1),prot(:,3),'r*',prot(:,1),prot(:,5),'b*');
axis([0 Nmax+10 0 1.1*max(prot(:,5))]); 
set(gca,'fontsize',14);
grid on;
xlabel('n');
ylabel('Laufzeit [s]');
legend('Broyden-Verfahren','Newton-Verfahren',4);

print -dpsc2 'bigbroydentiming.eps';
