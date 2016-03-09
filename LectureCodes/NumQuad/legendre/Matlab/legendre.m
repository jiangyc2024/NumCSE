function V= legendre(n,x)
V = ones(size(x)); V = [V; x];
for j=1:n-1
  V = [V; ((2*j+1)/(j+1)).*x.*V(end,:) - j/(j+1)*V(end-1,:)]; end
