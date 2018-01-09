%Script zu Aufgabe 1:

clear all
close all

t = (0.1:0.1:1)';
f = [100 34 17 12 9 6 5 4 4 2]';

gamma1 = data_fit_normal(t,f);
gamma2 = data_fit_qr(t,f);

% part c: create plot data:
A = make_A(t);
y1 = A*gamma1;
y2 = A*gamma2;

tl = (0.1:0.01:1)';
Al = make_A(tl);
yl1 = Al*gamma1;
yl2 = Al*gamma2;

figure(1), clf
subplot(1,2,1)
semilogy(tl, yl1, 'r', tl, yl2, 'b', t, f, 'k*')
legend('normal equation','qr fitting','data set')
xlabel('t')
ylabel('y')
subplot(1,2,2)
semilogy(t, (y1-f).^2, 'r*', t, (y2 - f).^2, 'bo')
legend('normal equation','qr fitting')
xlabel('t')
ylabel('fitting error')

% part d:
gamma1 - gamma2

fprintf('Residual of qr - residual of normal: %e\n',...
	norm(A*gamma2 - f) - norm(A*gamma1 - f))
fprintf('condition number of A:    %e\n',cond(A))
fprintf('condition number of A^tA: %e\n',cond(A'*A))