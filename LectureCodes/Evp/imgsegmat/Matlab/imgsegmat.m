function [A,D] = imgsegmat(P)
P = double(P); [n,m] = size(P);
spdata = zeros(4*n*m,1); spidx  = zeros(4*n*m,2);
k = 1;
for ni=1:n
  for mi=1:m
    mni = (mi-1)*n+ni;
    if (ni-1>0), spidx(k,:) = [mni,mni-1];
      spdata(k) = Sfun(P(ni,mi),P(ni-1,mi));
      k = k + 1;
    end
    if (ni+1<=n), spidx(k,:) = [mni,mni+1];
      spdata(k) = Sfun(P(ni,mi),P(ni+1,mi));
      k = k + 1;
    end
    if (mi-1>0), spidx(k,:) = [mni,mni-n];
      spdata(k) = Sfun(P(ni,mi),P(ni,mi-1));
      k = k + 1;
    end
    if (mi+1<=m), spidx(k,:) = [mni,mni+n];
      spdata(k) = Sfun(P(ni,mi),P(ni,mi+1));
      k = k + 1;
    end
  end
end
% Efficient initialization, see Sect.~\ref{sec:spml}, Ex.~\ref{ex:spinit}
W = sparse(spidx(1:k-1,1),spidx(1:k-1,2),spdata(1:k-1),n*m,n*m);
D = spdiags(full(sum(W')'),[0],n*m,n*m);
A = D-W;
 