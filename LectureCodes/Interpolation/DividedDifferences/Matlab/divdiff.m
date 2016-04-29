% 16.06.2009        divdiff.m
%
% Divided differences (two different implementations)
% t = node set (mutually different)
% y = nodal values
% Returns the coefficients of the polynomial in Newton basis

function y = divdiff(t,y)
n = length(y)-1;

% Non recursive implementation: backsubstitution
% for l=0:n-1
%   for j=l+1:n
%     y(j+1) = (y(j+1)-y(l+1))/(t(j+1)-t(l+1));
%   end
% end

% Recursive implementation
if (n > 0)
  y(1:n) = divdiff(t(1:n),y(1:n));
  for j=0:n-1
    y(n+1) = (y(n+1)-y(j+1))/(t(n+1)-t(j+1));
  end
end

