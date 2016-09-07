delimiterIn = ' ';
headerlinesIn = 0;
times = importdata('effmatmult.dat', delimiterIn, headerlinesIn);

% log-scale plot for investigation of asymptotic complexity
figure('name','dottenstiming');
loglog(times(:,1),times(:,3),'m+',...
       times(:,1),times(:,2),'r*',...
       times(:,1),times(:,1)*times(1,2)/times(1,1),'k-',...
       times(:,1),(times(:,1).^2)*times(2,2)/(times(2,1)^2),'b--');
xlabel('{\bf problem size n}','fontsize',14);
ylabel('{\bf average runtime (s)}','fontsize',14);
title('Timings for rank 1 matrix-vector multiplications');
legend('slow evaluation','efficient evaluation',...
       'O(n)','O(n^2)','location','northwest');
print -depsc2 './dottenstiming.eps';


