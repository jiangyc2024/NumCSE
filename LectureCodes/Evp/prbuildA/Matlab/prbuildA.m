function A = prbuildA(G,d)
N = size(G,1); 
l = full(sum(G)); idx = find(l>0);
s = zeros(N,1); s(idx) = 1./l(idx);
ds = ones(N,1)/N; ds(idx) = d/N;
A = ones(N,1)*ones(1,N)*diag(ds) + (1-d)*G*diag(s);
