delimiterIn = ' ';
headerlinesIn = 0;
times = importdata('arrowsystiming.dat', delimiterIn, headerlinesIn);

%% Dense plotting
figure('name','arrowsystiming_cpp');
loglog(times(:,1),times(:,2),'b-*',times(:,1),times(:,3),'r-+');
xlabel('{\bf matrix size n}','fontsize',14); 
ylabel('{\bf runtime [s]}','fontsize',14);
legend('arrowsys slow','arrowsys fast','location','northwest');

print -depsc2 './arrowsystiming_cpp.eps';

%% Sparse plotting
figure('name','arrowsystiming_sparse_cpp');
loglog(times(:,1),times(:,2),'b-*',times(:,1),times(:,3),'r-+',...
    times(:,1),times(:,4),'m-',times(:,1),times(:,5),'g-x');
xlabel('{\bf matrix size n}','fontsize',14); 
ylabel('{\bf runtime [s]}','fontsize',14);
legend('arrowsys slow','arrowsys fast','arrowsys SparseLU','arrowsys iterative ','location','northwest');

print -depsc2 './arrowsystiming_sparse_cpp.eps';