function comruleerrs()
% Numerical quadrature on [0,1]

figure('Name','1/(1+(5t)^2)');
exact = atan(5)/5;
trres = trapezoidal(inline('1./(1+(5*x).^2)'),0,1,1:200);
smres = simpson(inline('1./(1+(5*x).^2)'),0,1,1:200);
loglog(trres(:,1),abs(trres(:,2)-exact),'r+-',...
       smres(:,1),abs(smres(:,2)-exact),'b+-',...
       trres(:,1),trres(:,1).^2*(trres(1,2)/trres(1,1)^2),'r--',...
       smres(:,1),smres(:,1).^4*(smres(1,2)/smres(1,1)^2),'b--');
set(gca,'fontsize',12);
title('numerical quadrature of function 1/(1+(5t)^2)','fontsize',14);
xlabel('{\bf meshwidth}','fontsize',14);
ylabel('{\bf |quadrature error|}','fontsize',14);
legend('trapezoidal rule','Simpson rule','O(h^2)','O(h^4)',2);
axis([1/300 1 10^(-15) 1]);
trp1 = polyfit(log(trres(end-100:end,1)),log(abs(trres(end-100:end,2)-exact)),1)
smp1 = polyfit(log(smres(end-100:end,1)),log(abs(smres(end-100:end,2)-exact)),1)
print -dpsc2 '../PICTURES/compruleerr1.eps';

figure('Name','sqrt(t)');
exact = 2/3;
trres = trapezoidal(inline('sqrt(x)'),0,1,1:200);
smres = simpson(inline('sqrt(x)'),0,1,1:200);
loglog(trres(:,1),abs(trres(:,2)-exact),'r+-',...
       smres(:,1),abs(smres(:,2)-exact),'b+-',...
       trres(:,1),trres(:,1).^(1.5)*(trres(1,2)/trres(1,1)^2),'k--');
set(gca,'fontsize',14);
title('numerical quadrature of function sqrt(t)','fontsize',14);
xlabel('{\bf meshwidth}','fontsize',14);
ylabel('{\bf |quadrature error|}','fontsize',14);
legend('trapezoidal rule','Simpson rule','O(h^{1.5})',2);    
axis([1/300 1 10^(-7) 1]);
trp2 = polyfit(log(trres(end-100:end,1)),log(abs(trres(end-100:end,2)-exact)),1)
smp2 = polyfit(log(smres(end-100:end,1)),log(abs(smres(end-100:end,2)-exact)),1)
print -dpsc2 '../PICTURES/compruleerr2.eps';



