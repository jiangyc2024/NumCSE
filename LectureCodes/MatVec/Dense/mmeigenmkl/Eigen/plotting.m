delimiterIn = ' ';
headerlinesIn = 0;
times = importdata('eigen_sequential.dat', delimiterIn, headerlinesIn);
tmp1 = importdata('eigen_parallel.dat', delimiterIn, headerlinesIn);
tmp2 = importdata('mkl_sequential.dat', delimiterIn, headerlinesIn);
tmp3 = importdata('mkl_parallel.dat', delimiterIn, headerlinesIn);
times = [times, tmp1(:,2), tmp2(:,2), tmp3(:,2)];


figure('name','mmeigenmkl_1');
loglog(times(:,1),times(:,2),'r-+',times(:,1),times(:,3),'m-*',...
        times(:,1),times(:,4),'b-^', times(:,1),times(:,5),'k-o' );
xlabel('matrix size n','fontsize',14);
ylabel('execution time [s]','fontsize',14);
legend('Eigen sequential', 'Eigen parallel', 'MKL sequential','MLK parallel','location','northwest');

print -depsc2 './mmeigenmkl_1.eps';


figure('name','mmeigenmkl_2');
loglog(times(:,1),times(:,2)./times(:,1).^3,'r-+',times(:,1),times(:,3)./times(:,1).^3,'m-*',...
        times(:,1),times(:,4)./times(:,1).^3,'b-^', times(:,1),times(:,5)./times(:,1).^3,'k-o' );
xlabel('matrix size n','fontsize',14);
ylabel('execution time divided by n^3 [s]','fontsize',14);
legend('Eigen sequential', 'Eigen parallel', 'MKL sequential','MLK parallel','location','northwest');

print -depsc2 './mmeigenmkl_2.eps';

