function numquaderrs()
% Numerical quadrature on [0,1]
N = 20;

figure('Name','1/(1+(5t)^2)');
exact = atan(5)/5;
eqdres = numquad(inline('1./(1+(5*x).^2)'),0,1,N,'equidistant');
chbres = numquad(inline('1./(1+(5*x).^2)'),0,1,N,'Chebychev');
gaures = numquad(inline('1./(1+(5*x).^2)'),0,1,N,'Gauss');
semilogy(eqdres(:,1),abs(eqdres(:,2)-exact),'b+-',...
     chbres(:,1),abs(chbres(:,2)-exact),'m+-',...
     gaures(:,1),abs(gaures(:,2)-exact),'r+-');
set(gca,'fontsize',12);
title('Numerical quadrature of function  1/(1+(5t)^2)');
xlabel('{\bf Number of quadrature nodes}','fontsize',14);
ylabel('{\bf |quadrature error|}','fontsize',14);
legend('Equidistant Newton-Cotes quadrature',...
       'Clenshaw-Curtis quadrature',...
       'Gauss quadrature',3);
eqdp1 = polyfit(eqdres(:,1),log(abs(eqdres(:,2)-exact)),1)
chbp1 = polyfit(chbres(:,1),log(abs(chbres(:,2)-exact)),1)
gaup1 = polyfit(gaures(:,1),log(abs(gaures(:,2)-exact)),1)
print -depsc2 '../PICTURES/numquaderr1.eps';

figure('Name','sqrt(t)');
exact = 2/3;
eqdres = numquad(inline('sqrt(x)'),0,1,N,'equidistant');
chbres = numquad(inline('sqrt(x)'),0,1,N,'Chebychev');
gaures = numquad(inline('sqrt(x)'),0,1,N,'Gauss');
loglog(eqdres(:,1),abs(eqdres(:,2)-exact),'b+-',...
     chbres(:,1),abs(chbres(:,2)-exact),'m+-',...
     gaures(:,1),abs(gaures(:,2)-exact),'r+-');
set(gca,'fontsize',12);
axis([1 25 0.000001 1]);
title('Numerical quadrature of function sqrt(t)');    
xlabel('{\bf Number of quadrature nodes}','fontsize',14);
ylabel('{\bf |quadrature error|}','fontsize',14);
legend('Equidistant Newton-Cotes quadrature',...
       'Clenshaw-Curtis quadrature',...
       'Gauss quadrature',1);
eqdp2 = polyfit(log(eqdres(:,1)),log(abs(eqdres(:,2)-exact)),1)
chbp2 = polyfit(log(chbres(:,1)),log(abs(chbres(:,2)-exact)),1)
gaup2 = polyfit(log(gaures(:,1)),log(abs(gaures(:,2)-exact)),1)
print -depsc2 '../PICTURES/numquaderr2.eps';



