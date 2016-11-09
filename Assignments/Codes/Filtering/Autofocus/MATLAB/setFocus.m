function B = setFocus (f)
f0 = 2.0;  % true value of the focus (not known a-priori)
P =  double(imread('eth.pbm'));
ep = max(abs(f - f0),eps);
M = 16;
l = -M:M;
k = -M:M;
[L,K] = meshgrid(l,k);
S = 1 ./ (1 + (L.^2 + K.^2)/ep);
S = S / sum(sum(S));
B = conv2(P,S,'same');

