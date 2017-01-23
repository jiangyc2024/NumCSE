n = 1000;

% First setting

s = ones(n,1); % Exact solution
b = (1:n)';
y = b-s;
a = 1/sqrt(sum(y))*y;
A = (eye(n) + a*a');

% Second setting
% 
% s =  ones(n,1); % Exact solution
% A = gallery('minij',n);
% b = diag(s)*A*s;

% Non-linear function

F = @(x) (diag(x)*A*x - b);
DF = @(x) (diag(x)*A + diag(A*x));

% Initial guess
a = 2; b = 4; h = (b-a)/n;
x0 = (a:h:b-h)';
fprintf('Residual F(s) = %e\n',norm(F(s)));
fprintf('Initial Residual F(x) = %e\n',norm(F(x0)));

tol = 2*n*(1.0e-5);


% Convergence history
[x,res] = updbroyd(F,x0,DF(x0),tol,20,s);
[x,nres] = newton(x0,F,DF,size(res,1),s);
[x,sres] = simplenewt(x0,F,DF,size(res,1),s);

figure('name','Convergence of Broydens method');
semilogy(res(:,1),res(:,3),'m.',...
	 res(:,1),res(:,4),'r-*',...
	 nres(:,1),nres(:,3),'k.',...
	 nres(:,1),nres(:,4),'b-*',...
	 sres(:,1),sres(:,4),'c--*');
axis([0 size(res,1) 1.0e-13 1.0e5]);
set(gca,'fontsize',12);
grid on;
title(sprintf('n = %d, tol = %e',n,tol));
xlabel('Iterationsschritt');
ylabel('Normen');
legend('Broyden: ||F(x^{(k)})||','Broyden: Fehlernorm',...
       'Newton: ||F(x^{(k)})||','Newton: Fehlernorm',...
       'Newton (vereinfacht)',4);
print -dpsc2 'bigbroydencvg.eps';

