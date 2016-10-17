function plotspaicg ()
MAXIT = 1000;
ns = 2.^(1:4);
itCG = ones(1,length(ns));
itPCG = ones(1,length(ns));
Ns = zeros(1,length(ns));
for i = 1:length(ns)
    n = ns(i);
    A = kron( spdiags(exp((1:n*n)/(n*n))',0,n*n,n*n), ...
              gallery('tridiag',n*n,-1,2,-1) );
    Ns(i) = size(A,1);
    b = ones(Ns(i),1);
    [~,~,~,itCG(i)] = pcg(A,b,[],MAXIT);
    B = spai(A);
    [~,~,~,itPCG(i)] = pcg(A,b,[],MAXIT,@(x)(B+B')/2*x);
end
loglog(Ns, itCG, 'r-o', Ns, itPCG, 'b-x');
legend('plain CG', 'preconditioned CG', 'location', 'SouthEast');
title('{\bf Number of CG iterations vs. the system size N}');
xlabel('{\bf N}');
ylabel('{\bf CG iters}');
print -depsc '../PICTURES/plotspaicg.eps'
