spy(M);
[R,P] = lu(M); spy(R);
r = symrcm(M);
[R,P] = chol(M(r,r)); spy(R);
m = symamd(M);
[R,p] = chol(M(m,m)); spy(R);