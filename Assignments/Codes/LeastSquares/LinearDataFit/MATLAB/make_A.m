function A = make_A(t)


A(:,1) = 1./t;
A(:,2) = 1./t.^2;
A(:,3) = exp(-(t-1));
A(:,4) = exp(-2*(t-1));