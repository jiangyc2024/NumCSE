function C = blur(P,S)
[m,n] = size(P); [M,N] = size(S);
if (M ~= N), error('S not quadratic'); end
L = (M-1)/2; C = zeros(m,n);
for l=1:m, for j=1:n
    s = 0;
    for k=1:(2*L+1), for q=1:(2*L+1)
	kl = l+k-L-1;
	if (kl < 1), kl = kl + m; end
	if (kl > m), kl = kl - m; end
	jm = j+q-L-1;
	if (jm < 1), jm = jm + n; end
	if (jm > n), jm = jm - n; end
	s = s+P(kl,jm)*S(k,q);
      end, end
    C(l,j) = s;
  end, end
