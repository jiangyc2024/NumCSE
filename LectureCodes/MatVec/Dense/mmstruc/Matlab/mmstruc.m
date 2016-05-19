n = 100; A = [diag(1:n-1), (1:n-1)'; (1:n) ]; B = A(n:-1:1,:);
C = A*A; D = A*B;
figure; spy(A,'r'); axis off; print -depsc2 '../PICTURES/Aspy.eps';
figure; spy(B,'r'); axis off; print -depsc2 '../PICTURES/Bspy.eps';
figure; spy(C,'r'); axis off; print -depsc2 '../PICTURES/Cspy.eps';
figure; spy(D,'r'); axis off; print -depsc2 '../PICTURES/Dspy.eps';
