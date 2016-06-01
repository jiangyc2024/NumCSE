delimiterIn = ' ';
headerlinesIn = 0;
times = importdata('nomkl_gauss.dat', delimiterIn, headerlinesIn);
tmp1 = importdata('mkl_sequential.dat', delimiterIn, headerlinesIn);
tmp2 = importdata('mkl_parallel.dat', delimiterIn, headerlinesIn);
times = [times, tmp1(:,2), tmp2(:,2)];

figure('name','gausstiming');
loglog(times(:,1),times(:,2),'r-+',times(:,1),times(:,3),'m-*',...
        times(:,1),times(:,4),'b-^', times(:,1),times(:,5),'g-<',...
       times(:,1),times(:,1).^3*(1E-6/(times(1,1)^3)),'k-');
xlabel('matrix size n','fontsize',14);
ylabel('execution time [s]','fontsize',14);
legend('Eigen lu() solver', 'gausselimsolve','MLK solver sequential','MLK solver parallel','O(n^3)','location','northwest');

print -depsc2 './gausstiming.eps';
