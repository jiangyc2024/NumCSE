function z = sdirkStep(z0,h,gamma)
% SDIRKSTEPVER1 calculate one SDIRK step for the differential equation
% in exercise sheet 9, task 5.
%
% See also SDIRK, SDIRKCONV, SDIRKSTEPVER2.

  A = [ 0  1
       -1 -1];
  
  %Calculate k1 and k2
  k = nan(2,size(z0,1));

  k(1,:) = ( eye(size(A)) - h*gamma*A )\(A*z0);
  k(2,:) = ( eye(size(A)) - h*gamma*A )\(A*z0 + h*(1-2*gamma)*A*k(1,:).');
  
  
  % perform one SDIRK-step
  z = z0 + h/2*( k(1,:).' + k(2,:).' );