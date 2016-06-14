% Curing Wilkinson's counterexample by random perturbation
% Theory: Spielman and Teng
res = [];
for n=10:10:200
  % Build Wilkinson matrix
  A = [tril(-ones(n,n-1))+2*[eye(n-1);
       zeros(1,n-1)],ones(n,1)];
  % imposed solution
  x = ((-1).^(1:n))';
  relerr = norm(A\(A*x)-x)/norm(x);
  % Randomly perturbed Wilkinson matrix by matrix with iid
  % $N(0,{\mathrm{eps}})$ distributed entries
  Ap = A + eps*randn(size(A));
  relerrp = norm(Ap\(A*x)-x)/norm(x);
  res = [res; n relerr relerrp];
end
semilogy(res(:,1),res(:,2),'m-*',res(:,1),res(:,3),'r-+');
xlabel('matrix size n','fontsize',14);
ylabel('relative error','fontsize',14);
legend('unperturbed matrix','randn perturbed matrix','location','west');

print -depsc2 '../PICTURES/wilkpert.eps';
