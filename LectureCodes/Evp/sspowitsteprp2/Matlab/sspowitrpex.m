function sspowitrpex(d,maxit)
% monitor power iteration with Ritz projection for computing
% the two largest (in modulus) eigenvalues and associated eigenvectors
% of a symmetric matrix with prescribed eigenvalues passed in \texttt{d}
if (nargin < 10), maxit = 20; end
if (nargin < 1), d = (1:10)'; end
% Generate matrix
n = length(d);
Z = diag(sqrt(1:n),0) + ones(n,n);
[Q,R] = qr(Z);   % generate orthogonal matrix
A = Q*diag(d,0)*Q'; % \Blue{$\VA=\VA^T$} with spectrum \Blue{$\sigma(\VA) =\{d_1,\ldots,d_n\}$}
% Compute ``exact'' eigenvectors and eigenvalues
[V,D] = eig(A); [d,idx] = sort(diag(D)),
v_ex = V(:,idx(n)); w_ex = V(:,idx(n-1)); 
lv_ex = d(n); lw_ex = d(n-1);

v = ones(n,1); w = (-1).^v; % (Arbitrary) initial guess for eigenvectors
v = v/norm(v); w = w/norm(w);
result = [];
A
for k=1:maxit
  v_new = A*v; w_new = A*w; % ``power iteration'', \emph{cf.} \eqref{eq:dirpotmeth}
  [Q,R] = qr([v_new,w_new],0); % orthogonalization, \emph{cf.} Rem.~\ref{rem:QRorth}
  Q
  break
  [U,D] = eig(Q'*A*Q); % Solve Ritz projected eigenvalue problem
  [ev,idx] = sort(abs(diag(D))), % Sort eigenvalues
  w = Q*U(:,idx(1)); v = Q*U(:,idx(2)); % Recover approximate eigenvectors
  
  % Record errors in eigenvalue and eigenvector approximations. Note that the 
  % direction of the eigenvectors is not specified.
  result = [result; k, abs(ev(2)-lv_ex), abs(ev(1)-lw_ex), ...
            min(norm(v-v_ex),norm(v+v_ex)), min(norm(w-w_ex),norm(w+w_ex))];
end

figure('name','sspowitrp');
semilogy(result(:,1),result(:,2),'m-+',...
     result(:,1),result(:,3),'r-*',...
     result(:,1),result(:,4),'k-^',...
     result(:,1),result(:,5),'b-p');
title('d = [0.5*(1:8),9.5,10]');
xlabel('{\bf power iteration step}','fontsize',14);
ylabel('{\bf error}','fontsize',14);
legend('error in \lambda_n','error in \lambda_n-1','error in v','error in w','location','northeast');
print -depsc2 '../PICTURES/sspowitrpcvg.eps';

rates = result(2:end,2:end)./result(1:end-1,2:end);
figure('name','rates(sspowitrp)');
plot(result(2:end,1),rates(:,1),'m-+',...
     result(2:end,1),rates(:,2),'r-*',...
     result(2:end,1),rates(:,3),'k-^',...
     result(2:end,1),rates(:,4),'b-p');
axis([0 maxit 0 1]);
title('d = [0.5*(1:8),9.5,10]');
xlabel('{\bf power iteration step}','fontsize',14);
ylabel('{\bf error quotient}','fontsize',14);
legend('error in \lambda_n','error in \lambda_n-1','error in v','error in w','location','northeast');
print -depsc2 '../PICTURES/sspowitrpcvgrates.eps';


     
