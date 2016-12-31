n = 1000;
A = n*eye(n) + rand(n,n);
b = rand(n,1);

tic;
for k=1:100
  x = A\b;
end
t = toc;
fprintf('time (backslash) = %e\n',t);

tic;
[L,U] = lu(A);
for k=1:100
  y = U\(L\b);
end
t = toc;
fprintf('time (lu) = %e\n',t);

tic;
Ai = inv(A);
for k=1:100
  y = Ai*b;
end
t = toc;
fprintf('time (invmul) = %e\n',t);
