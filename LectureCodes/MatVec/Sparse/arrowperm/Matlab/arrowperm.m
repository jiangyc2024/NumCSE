n = 10; A = [ n+1, (n:-1:1) ; ones(n,1), eye(n,n)];
% Permutation matrix ($\to$ Def.~\ref{def:permmat}) encoding cyclic permutation
P = [zeros(n,1),eye(n);1,zeros(1,n)];

figure('name','A');    
spy(A,'r.'); print -depsc  '../PICTURES/InvArrowSpy.eps';
figure('name','PAPT');
spy(P*A*P','r.'); print -depsc '../PICTURES/ArrowSpy.eps';

