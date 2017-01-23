function [x,prot] = updbroyd(F,x,J,tol,maxit,sol)
% Good Broyden rank-1-update quasi-Newton method
% Straightforward implementation for small problems
% F = Non-linear mapping in n dimensions
% x = initial guess
% J = initial guess for Jacobi matrix at x0
% tol = tolerance for termination

if (nargin < 6), sol = []; end
if (nargin < 5), maxit = 20; end

n = length(x);
% Initial step 
Ji = inv(J);
% Storing ineverses of approximate Jacobians
Jstor(:,:,1) = Ji; x0 = x;

k = 1; tol = tol^2; s = Ji*F(x); 
x = x - s; f = F(x); sn = dot(s,s);
fprintf('Iteration(UPD) %d: |s| = %e, |F(x)| = %e \n',k,sqrt(sn),norm(f));

% Storing update vectors
Sstor = [s;sn];

% Keeping a record of convergence history
prot = [k,sqrt(sn),norm(f)];

if (~isempty(sol)), prot = [prot, norm(x-sol)]; end
while ((sn > tol) && (k < 10))
  Ju = inv(J);
  w = Ju*f;
  for l=2:k
    w = w+Sstor(1:end-1,l)*(Sstor(1:end-1,l-1)'*w)/Sstor(end,l-1);
    Ju = (eye(n) + (Sstor(1:end-1,l)*Sstor(1:end-1,l-1)')/Sstor(end,l-1))*Ju;  
  end
  v = Ji*f;
  fprintf('|Ji-Ju| = %e, |w-v| = %e\n',norm(Ji-Ju,'fro'),norm(w-v));
  Ji = (eye(n)+v*s'/(sn-s'*v))*Ji;
  s = Ji*f; x = x - s; 
  f = F(x); sn = dot(s,s);
  k = k+1;
  fprintf('Iteration(UPD) %d: |s| = %e, |F(x)| = %e \n',k,sqrt(sn),norm(f));
  Jstor(:,:,k) = Ji;
  Sstor = [Sstor,[s;sn]];
  res = [k,sqrt(sn),norm(f)];
  if (~isempty(sol)), res = [res, norm(x-sol)]; end
  prot = [prot; res];
end

% Redoing iteration 
f = F(x0);
for l=1:k
  s = Jstor(:,:,l)*f;
  x0 = x0-s;
  f = F(x0);
  fprintf('Retrace(UPD) %d: |s| = %e, |F(x)| = %e \n',l,norm(s),norm(f));
end

% Checking update of Jacobians
Ju = Jstor(:,:,1);
for l=2:k
  Ju = (eye(n) + (Sstor(1:end-1,l)*Sstor(1:end-1,l-1)')/Sstor(end,l-1))*Ju;
  fprintf('Step %d: diff = %e\n',l,norm(Jstor(:,:,l)-Ju,'fro'));
end
