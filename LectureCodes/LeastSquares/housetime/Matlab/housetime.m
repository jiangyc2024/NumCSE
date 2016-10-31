% Timing QR factorizations

K = 3; r = [];
for n=2.^(2:6)
  m = n*n;
  
  A = (1:m)'*(1:n) + [eye(n);ones(m-n,n)];
  t1 = 1000; for k=1:K, tic; [Q,R] = qr(A); t1 = min(t1,toc); clear Q,R; end
  t2 = 1000; for k=1:K, tic; [Q,R] = qr(A,0); t2 = min(t2,toc); clear Q,R; end
  t3 = 1000; for k=1:K, tic; R = qr(A); t3 = min(t3,toc); clear R; end
  r = [r; n , m , t1 , t2 , t3];
end


