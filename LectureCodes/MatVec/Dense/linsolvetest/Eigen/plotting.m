delimiterIn = ' ';
headerlinesIn = 0;
times = importdata('linsolvetest.dat', delimiterIn, headerlinesIn);

figure('name','linsolvetest_cpp');
loglog(times(:,1),times(:,2),'r+',times(:,1),times(:,3),'m*');
xlabel('{\bf matrix size n}','fontsize',14); 
ylabel('{\bf runtime for direct solver [s]}','fontsize',14);
legend('naive lu()','triangularView lu()','location','northwest');

print -depsc2 './linsolvetest_cpp.eps';
