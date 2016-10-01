function [G,x] = planerot(a)
% plane Givens rotation avoiding cancellation.
if (a(2) ~= 0)
  r = norm(a);  G = [a'; -a(2) a(1)]/r; x = [r; 0];
else G = eye(2);
end
