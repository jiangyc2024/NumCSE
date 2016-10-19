function  y = tridiagqr(c,d,e,b)
n = length(d); t = norm(d)+norm(e)+norm(c);
for k=1:n-1
  [R,z] = planerot([d(k);e(k)]);
  if (abs(z(1))/t < eps), error('Matrix singular'); end;
  d(k) = z(1); b(k:k+1) = R*b(k:k+1);
  Z = R*[c(k), 0;d(k+1), c(k+1)];
  c(k) = Z(1,1); d(k+1) = Z(2,1);
  e(k) = Z(1,2); c(k+1) = Z(2,2);
end
% Backsubstitution; do directly on vectors
A = spdiags([d,[0;c(1:end-1)],[0;0;e(1:end-2)]],[0 1 2],n,n);
y = A\b;
