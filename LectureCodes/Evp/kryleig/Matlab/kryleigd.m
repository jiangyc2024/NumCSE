% MATLAB script for ex:kryleig
n=100; M=gallery('tridiag',-0.5*ones(n-1,1),2*ones(n,1),-1.5*ones(n-1,1));
[Q,R]=qr(M); A=Q'*diag(1:n)*Q; % eigenvalues \Blue{$1,2,3,\ldots,100$}

[V0,D0] = eig(A); d0 = sort(diag(D0)); 
ev0 = [d0(1:3);d0(end-2:end)];

errev = [];
for m=6:30
  [V,D] = kryleig(A,m); d = sort(diag(D));
  ev =  [d(1:3);d(end-2:end)];
  errev = [errev; m , ev', abs((ev-ev0)')];
end

figure; 
plot(errev(:,1),errev(:,2),'r+',errev(:,1),errev(:,3),'m+',errev(:,1),errev(:,4),'b+');
axis([5 31 -1 41]);
xlabel('{\bf dimension m of Krylov space}','fontsize',14);
ylabel('{\bf Ritz value}','fontsize',14);
legend('\mu_1','\mu_2','\mu_3','location','northeast');

print -depsc2 '../PICTURES/kryleigmin.eps';

figure;
plot(errev(:,1),errev(:,7),'r^',errev(:,1),errev(:,6),'m^',errev(:,1),errev(:,5),'b^');
axis([5 31 70 101]);
xlabel('{\bf dimension m of Krylov space}','fontsize',14);
ylabel('{\bf Ritz value}','fontsize',14);
legend('\mu_{m}','\mu_{m-1}','\mu_{m-2}','location','southeast');

print -depsc2 '../PICTURES/kryleigmax.eps';


figure;
semilogy(errev(:,1),errev(:,8),'r+',errev(:,1),errev(:,9),'m+',errev(:,1),errev(:,10),'b+');
ax = axis; axis([5 31 ax(3) ax(4)]);
xlabel('{\bf dimension m of Krylov space}','fontsize',14);
ylabel('{\bf error of Ritz value}','fontsize',14);
legend('|\lambda_1-\mu_1|','|\lambda_2-\mu_2|','|\lambda_2-\mu_3|',...
       'location','northeast');

print -depsc2 '../PICTURES/kryleigminerr.eps';

figure;
semilogy(errev(:,1),errev(:,13),'r^',errev(:,1),errev(:,12),'m^',errev(:,1),errev(:,11),'b^');
ax = axis; axis([5 31 ax(3) ax(4)]);
xlabel('{\bf dimension m of Krylov space}','fontsize',14);
ylabel('{\bf Ritz value}','fontsize',14);
legend('|\lambda_m-\mu_{m}|','|\lambda_{m-1}-\mu_{m-1}|','|\lambda_{m-2}-\mu_{m-2}|','location','northeast');

print -depsc2 '../PICTURES/kryleigmaxerr.eps';




