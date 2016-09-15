function ModNewtOrder
a = 0.123 ;
F = @(x) atan(x) - a ;
DF = @(x) 1./(1 + x.^2) ;
x_exact = tan(a);
% x_exact = fzero(F,x0);
x0 = 5;

x = x0;
err = abs(x0 - x_exact);
it = 0;
while (err(end) > eps   &&   it < 100)
    x_new = ModNewtStep (x(end), F, DF);
    x = [x, x_new];
    err = [ err, abs(x_new - x_exact) ];
    it = it+1;
end
% logarithm of the error and ratios to estimate order of conv.:
log_err = log(err);
emp_orders = (log_err(3:end) - log_err(2:end-1)) ./ ...
    (log_err(2:end-1) - log_err(1:end-2));

close all; figure;
semilogy(1:length(err), err,'o-');
xlabel('{\bf Iteration number}','fontsize',14);
ylabel('{\bf Error}','fontsize',14);
print -depsc2  '../PICTURES/ModNewtOrder.eps';

disp('Errors and empirical orders of conv.:')
[err', [emp_orders';0;0]]
