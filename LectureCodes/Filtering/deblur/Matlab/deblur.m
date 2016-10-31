function D = deblur(C,S,tol)
[m,n] = size(C); [M,N] = size(S);
if (M ~= N), error('S not quadratic'); end
L = (M-1)/2; Spad = zeros(m,n);
% Zero padding 
Spad(1:L+1,1:L+1) = S(L+1:end,L+1:end);
Spad(m-L+1:m,n-L+1:n) = S(1:L,1:L);
Spad(1:L+1,n-L+1:n) = S(L+1:end,1:L);
Spad(m-L+1:m,1:L+1) = S(1:L,L+1:end);
% Inverse of blurring operator
SF = fft2(Spad); 
% Test for invertibility
if (nargin < 3), tol = 1E-3; end
if (min(min(abs(SF))) < tol*max(max(abs(SF))))
  error('Deblurring impossible'); 
end
% DFT based deblurring
D = fft2(ifft2(C)./SF);
