function [V,res] = lanczosev(A,k,j,w)
V = [w/norm(w)];
alpha = zeros(k,1);
beta = zeros(k,1);
res = [];
for l=1:k
  vt = A*V(:,l);
  if (l>1)
    vt = vt - beta(l-1)*V(:,l-1);
  end
  alpha(l) = dot(V(:,l),vt);
  vt = vt - alpha(l)*V(:,l);
  beta(l) = norm(vt);
  V = [V, vt/beta(l)];
  T2 = spdiags([[beta(1:l-1);0],alpha(1:l),[0;beta(1:l-1)]],[-1,0,1],l,l);
  ev = sort(eig(full(T2)));
  if (j > l)
    ev = [zeros(j-l,1); ev]; 
    res = [res; l ev']; 
  else 
    res = [res; l ev(end-j+1:end)'];
  end
end
