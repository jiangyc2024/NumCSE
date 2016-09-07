delimiterIn = ' ';
headerlinesIn = 0;
timings = importdata('scaletiming.dat', delimiterIn, headerlinesIn);

figure('name','scaletimings');
loglog(timings(:,1),timings(:,2),'r*',...
       timings(:,1),timings(:,3),'b+');
xlabel('{\bf vector length n}','fontsize',14);
ylabel('{\bf time [s]}','fontsize',14);
title('Timings for different ways to do scaling');
legend('D.diagonal().cwiseProduct(x)','D*x','location','best');

print -depsc2 './scaletiming_m.eps';



