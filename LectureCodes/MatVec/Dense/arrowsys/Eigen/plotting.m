delimiterIn = ' ';
headerlinesIn = 0;
times = importdata('arrowsystiming.dat', delimiterIn, headerlinesIn);

figure('name','arrowsystiming_cpp');
loglog(times(:,1),times(:,2),'b-*',times(:,1),times(:,3),'r-+');
xlabel('{\bf matrix size n}','fontsize',14); 
ylabel('{\bf runtime [s]}','fontsize',14);
legend('arrowsys slow','arrowsys fast','location','northwest');

print -depsc2 './arrowsystiming_cpp.eps';
