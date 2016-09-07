delimiterIn = ' ';
headerlinesIn = 0;
times = importdata('mmtiming.dat', delimiterIn, headerlinesIn);

figure('name','mmtiming');
loglog(times(:,1),times(:,2),'r+-',...
     times(:,1),times(:,3),'m*-',...
     times(:,1),times(:,4),'b^-',...
     times(:,1),times(:,5),'kp-');
title('Timings: Different implementations of matrix multiplication');
xlabel('matrix size n','fontsize',14);
ylabel('time [s]','fontsize',14);
legend('loop implementation','dot-product implementation',...
       'matrix-vector implementation','Eigen matrix product',...
       'location','northwest');
print -depsc2 './mmtiming_cpp.eps';


