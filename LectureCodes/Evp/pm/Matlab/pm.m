% Demonstration of direct power method for Ex.~\ref{ex:pm}
maxit = 30; d = [1:10]'; n = length(d);
% Initialize the matrix \Blue{$\VA$}
S = triu(diag(n:-1:1,0)+ones(n,n));
A = S*diag(d,0)*inv(S);
% This calculates the exact eigenvalues (for error calculation)
[V,D] = eig(A);

k = find(d == max(abs(d)));
if (length(k) > 1), error('No single largest EV'); end;
ev = X(:,k(1)); ev = ev/norm(ev); ev
ew = d(k(1)); ew
if (ew < 0), sgn = -1; else sgn = 1; end

z = rand(n,1); z = z/norm(z);
s = 1;
res = []; 

% Actual direct power iteration
for i=1:maxit
  w = A*z; l = norm(w); rq = real(dot(w,z)); z = w/l;
  res = [res;i,l,norm(s*z-ev),abs(l-abs(ew)),abs(sgn*rq-ew)];
  s = s * sgn;
end

% Plot the result
semilogy(res(:,1),res(:,3),'r-*',res(:,1),res(:,4),'k-+',res(:,1),res(:,5),'m-o');
xlabel('{\bf iteration step k}','FontSize',14);
ylabel('{\bf errors}','FontSize',14);
print -deps2c '../PICTURES/pm1.eps';
