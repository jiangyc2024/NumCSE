% Driver for invstab

invstab;

figure('name','invtab');
loglog(result(:,1),result(:,2),'r+-',...
        result(:,1),result(:,3),'b-*',...
        result(:,1),result(:,4),'m^');
xlabel('{\bf \epsilon}','fontsize',14);
ylabel('{\bf relative residual}','fontsize',14);
legend('Gaussian elimination','multiplication with inverse',...
       'residual for inverse');

print -depsc2 '../PICTURES/invstab.eps';




