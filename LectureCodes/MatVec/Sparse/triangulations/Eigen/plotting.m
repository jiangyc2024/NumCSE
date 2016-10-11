delimiterIn = ' ';
headerlinesIn = 0;
times = importdata('timing.dat', delimiterIn, headerlinesIn);

exp = log(times(end-2,2)/times(end,2))/log(times(end-2,1)/times(end,1));

figure('name','Graph Laplacian');
loglog(times(:,1),times(:,2),'r-+',times(:,1),times(:,3),'m-*',...
        times(:,1),times(:,4),'b-^', times(:,1),times(:,5),'g-<',...
       times(:,1),times(:,6),'c-o', times(:,1), times(:,1).^exp / times(end,1).^exp * times(end,2),'k-');
xlabel('matrix size n','fontsize',14);
xlabel(' size of matrix A_{int}','fontsize',14);
ylabel('{\bf solution time [s]}','fontsize',14);
legend('Eigen SparseLU', 'Eigen SimplicialLDLT','Eigen ConjugateGradient','MKL PardisoLU','MKL PardisoLDLT',sprintf('O(n^{%.1f})', exp),'location','northwest');
print -depsc 'graphlaplsolvetiming_cpp.eps';

