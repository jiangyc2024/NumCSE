function [x,prot] = updbroyd(F,x,J,tol,maxit,sol)
% Good Broyden rank-1-update quasi-Newton method
% Straightforward implementation for small problems
% F = Non-linear mapping in n dimensions
% x = initial guess
% J = initial guess for Jacobi matrix at x0
% tol = tolerance for termination

if (nargin < 6), sol = []; end
if (nargin < 5), maxit = 20; end

k = 1; tol = tol^2; s = J\F(x); 
x = x - s; f = F(x); sn = dot(s,s);
fprintf('Iteration(UPD) %d: |s| = %e, |F(x)| = %e \n',k,sqrt(sn),norm(f));

% Storing update vectors
dx = [s];
dxn = [sn];

% Keeping a record of convergence history
prot = [k,sqrt(sn),norm(f)];

if (~isempty(sol)), prot = [prot, norm(x-sol)]; end
prot = [prot,1];
while ((sn > tol) && (k < maxit))
  w = J\f;
  for l=2:k
    w = w+dx(:,l)*(dx(:,l-1)'*w)/dxn(l-1);
  end
  z = s'*w; s = (1+z/(sn-z))*w; sn = s'*s; dx = [dx,s]; dxn = [dxn,sn];
  x = x - s; f = F(x); 
  k = k+1;
  fprintf('Iteration(UPD) %d: |s| = %e, |F(x)| = %e, theta = %f\n',...
	  k,norm(s),norm(f),norm(w)/sqrt(dxn(k-1)));
  res = [k,sqrt(sn),norm(f)];
  if (~isempty(sol)), res = [res, norm(x-sol)]; end
  prot = [prot; [res,norm(w)/sqrt(dxn(k-1))]];
end

%OUTPUT of debugging version:
% >> broydenmain
% Iteration 1: |s| = 1.285667, |F(x)| = 1.850483 
% Iteration 2: |s| = 1.204999, |F(x)| = 0.218124 
% Iteration 3: |s| = 0.059908, |F(x)| = 0.090416 
% Iteration 4: |s| = 0.019061, |F(x)| = 0.005713 
% Iteration 5: |s| = 0.007413, |F(x)| = 0.004380 
% Iteration 6: |s| = 0.009159, |F(x)| = 0.000370 
% Iteration 7: |s| = 0.000833, |F(x)| = 0.000006 
% Iteration 8: |s| = 0.000005, |F(x)| = 0.000000 
% Iteration 9: |s| = 0.000000, |F(x)| = 0.000000 
% Iteration(UPD) 1: |s| = 4.056678e-01, |F(x)| = 6.068983e-01 
% |Ji-Ju| = 0.000000e+00, |w-v| = 0.000000e+00
% Iteration(UPD) 2: |s| = 1.285667e+00, |F(x)| = 1.850483e+00 
% |Ji-Ju| = 5.102197e-15, |w-v| = 5.024296e-15
% Iteration(UPD) 3: |s| = 1.204999e+00, |F(x)| = 2.181238e-01 
% |Ji-Ju| = 1.553043e-14, |w-v| = 7.757919e-17
% Iteration(UPD) 4: |s| = 5.990776e-02, |F(x)| = 9.041630e-02 
% |Ji-Ju| = 1.567785e-14, |w-v| = 2.390459e-16
% Iteration(UPD) 5: |s| = 1.906113e-02, |F(x)| = 5.712561e-03 
% |Ji-Ju| = 1.151616e-14, |w-v| = 6.463446e-17
% Iteration(UPD) 6: |s| = 7.412905e-03, |F(x)| = 4.380372e-03 
% |Ji-Ju| = 1.383088e-14, |w-v| = 4.399083e-17
% Iteration(UPD) 7: |s| = 9.159427e-03, |F(x)| = 3.704975e-04 
% |Ji-Ju| = 2.318533e-14, |w-v| = 7.601796e-18
% Iteration(UPD) 8: |s| = 8.327204e-04, |F(x)| = 5.563950e-06 
% |Ji-Ju| = 2.168803e-14, |w-v| = 4.242781e-20
% Iteration(UPD) 9: |s| = 4.816393e-06, |F(x)| = 2.748540e-07 
% |Ji-Ju| = 2.184912e-14, |w-v| = 3.537951e-21
% Iteration(UPD) 10: |s| = 4.504276e-07, |F(x)| = 2.021961e-08 
% Retrace(UPD) 1: |s| = 4.056678e-01, |F(x)| = 6.068983e-01 
% Retrace(UPD) 2: |s| = 1.285667e+00, |F(x)| = 1.850483e+00 
% Retrace(UPD) 3: |s| = 1.204999e+00, |F(x)| = 2.181238e-01 
% Retrace(UPD) 4: |s| = 5.990776e-02, |F(x)| = 9.041630e-02 
% Retrace(UPD) 5: |s| = 1.906113e-02, |F(x)| = 5.712561e-03 
% Retrace(UPD) 6: |s| = 7.412905e-03, |F(x)| = 4.380372e-03 
% Retrace(UPD) 7: |s| = 9.159427e-03, |F(x)| = 3.704975e-04 
% Retrace(UPD) 8: |s| = 8.327204e-04, |F(x)| = 5.563950e-06 
% Retrace(UPD) 9: |s| = 4.816393e-06, |F(x)| = 2.748540e-07 
% Retrace(UPD) 10: |s| = 4.504276e-07, |F(x)| = 2.021961e-08 
% Step 2: diff = 5.102197e-15
% Step 3: diff = 1.553043e-14
% Step 4: diff = 1.567785e-14
% Step 5: diff = 1.151616e-14
% Step 6: diff = 1.383088e-14
% Step 7: diff = 2.318533e-14
% Step 8: diff = 2.168803e-14
% Step 9: diff = 2.184912e-14
% Step 10: diff = 2.005705e-14
