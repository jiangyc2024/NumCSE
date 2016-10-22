function sppowitdriver(d,maxit)
% monitor power iteration with orthogonal projection for finding
% the two largest (in modulus) eigenvalues and associated eigenvectors
% of a symmetric matrix with prescribed eigenvalues passed in \texttt{d}
if (nargin < 10), maxit = 20; end
if (nargin < 1), d = (1:10)'; end
% Generate matrix
n = length(d);
Z = diag(sqrt(1:n),0) + ones(n,n);
[Q,R] = qr(Z);   % generate orthogonal matrix
A = Q*diag(d,0)*Q'; % ``synthetic'' \Blue{$\VA=\VA^T$} with spectrum \Blue{$\sigma(\VA) =\{d_1,\ldots,d_n\}$}
% Compute ``exact'' eigenvectors and eigenvalues
[V,D] = eig(A); [d,idx] = sort(diag(D)),
v_ex = V(:,idx(n)); w_ex = V(:,idx(n-1)); 
lv_ex = d(n); lw_ex = d(n-1);

v = ones(n,1); w = (-1).^v; % (Arbitrary) initial guess for eigenvectors
v = v/norm(v); w = w/norm(w);
result = [];
for k=1:maxit
  v_new = A*v; w_new = A*w; % ``power iteration'', \emph{cf.} \eqref{eq:dirpotmeth}
  % Rayleigh quotients provide approximate eigenvalues
  lv = dot(v_new,v); lw = dot(w_new,w); 
  % orthogonalization, \emph{cf.} Gram-Schmidt orthogonalization \eqref{GRS}: \Blue{$\Vw\perp \Vv$}
  v = v_new/norm(v_new); w = w_new - dot(v,w_new)*v; w = w/norm(w);
  % Record errors in eigenvalue and eigenvector approximations. Note that the 
  % direction of the eigenvectors is not specified.
  result = [result; k, abs(lv-lv_ex), abs(lw-lw_ex), ...
            min(norm(v-v_ex),norm(v+v_ex)), min(norm(w-w_ex),norm(w+w_ex))];
end

figure('name','sspowit');
semilogy(result(:,1),result(:,2),'m-+',...
     result(:,1),result(:,3),'r-*',...
     result(:,1),result(:,4),'k-^',...
     result(:,1),result(:,5),'b-p');
title('d = [0.5*(1:8),9.5,10]');
xlabel('{\bf power iteration step}','fontsize',14);
ylabel('{\bf error}','fontsize',14);
legend('error in \lambda_n','error in \lambda_n-1','error in v','error in w','location','northeast');
print -depsc2 '../PICTURES/sspowitcvg1.eps';

rates = result(2:end,2:end)./result(1:end-1,2:end);
figure('name','rates');
plot(result(2:end,1),rates(:,1),'m-+',...
     result(2:end,1),rates(:,2),'r-*',...
     result(2:end,1),rates(:,3),'k-^',...
     result(2:end,1),rates(:,4),'b-p');
axis([0 maxit 0.5 1]);
title('d = [0.5*(1:8),9.5,10]');
xlabel('{\bf power iteration step}','fontsize',14);
ylabel('{\bf error quotient}','fontsize',14);
legend('error in \lambda_n','error in \lambda_n-1','error in v','error in w','location','southeast');
print -depsc2 '../PICTURES/sspowitcvgrates1.eps';


     
