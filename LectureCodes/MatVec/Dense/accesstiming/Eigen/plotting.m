delimiterIn = ' ';
headerlinesIn = 0;
times = importdata('timing_matlab.dat', delimiterIn, headerlinesIn);
tmp1 = importdata('timing_eigen.dat', delimiterIn, headerlinesIn);
times = [times, tmp1(:,2), tmp1(:,3)];

figure('name','accesstiming');
loglog(times(:,1),times(:,2),'r+', times(:,1),times(:,3),'m*', ...
    times(:,1),times(:,4),'bs', times(:,1),times(:,5),'k^');
xlabel('{\bf n}','fontsize',14); 
ylabel('{\bf runtime [s]}','fontsize',14); 
legend('A(:,j+1) = A(:,j+1) - A(:,j)','A(i+1,:) = A(i+1,:) - A(i,:)',...
        'eigen row access', 'eigen column access', ...
       'location','northwest');
print -depsc2 './accesstiming.eps';
