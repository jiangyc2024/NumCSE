% MATLAB script for timing different implementations of matrix multiplications
nruns = 3; times = [];
for n=2.^(2:10) % n = size of matrices
  fprintf('matrix size n = %d\n',n);
  A = rand(n,n); B = rand(n,n); C = zeros(n,n); 
  t1 = realmax;
  % loop based implementation (no BLAS)
  for l=1:nruns
    tic;
    for i=1:n, for j=1:n
        for k=1:n, C(i,j) = C(i,j) + A(i,k)*B(k,j); end
      end, end
    t1 = min(t1,toc);
    end
  t2 = realmax;
  % dot product based implementation (BLAS level 1)
  for l=1:nruns
    tic;
    for i=1:n
      for j=1:n, C(i,j) = dot(A(i,:),B(:,j));  end
    end
    t2 = min(t2,toc);   
  end
  t3 = realmax;
  % matrix-vector based implementation (BLAS level 2)
  for l=1:nruns
    tic;
    for j=1:n, C(:,j) = A*B(:,j); end
    t3 = min(t3,toc);
  end
  t4 = realmax;
  % BLAS level 3 matrix multiplication
  for l=1:nruns
    tic; C = A*B; t4 = min(t4,toc);
  end
  times = [ times; n t1 t2 t3 t4];
end
  
figure('name','mmtiming');
loglog(times(:,1),times(:,2),'r+-',...
     times(:,1),times(:,3),'m*-',...
     times(:,1),times(:,4),'b^-',...
     times(:,1),times(:,5),'kp-');
title('Timings: Different implementations of matrix multiplication');
xlabel('matrix size n','fontsize',14);
ylabel('time [s]','fontsize',14);
legend('loop implementation','dot product implementation',...
       'matrix-vector implementation','BLAS gemm (MATLAB *)',...
       'location','northwest');

print -depsc2 '../PICTURES/mvtiming.eps';