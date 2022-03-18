function x = GNradioactive
close all; tic;
% get measured data:   [t, m] column vectors
load decay.mat;

format shortE;

% CODE TEST: create perturbed vector  m  with known parameters
% xtest = [1.7,  0.4,  0.7,  0.2]';
% m = F(xtest, t, 0) + 0.01*randn(size(t));

% initial guess:
xi = [1,  1,  1, 0.1]';
x = xi;
GNupdate = 1;
GNiter = 0;
while GNupdate > 1e-14
    GNiter = GNiter+1;
    Fx = F(x, t, m);
    DFx = DF(x, t);
    % solve linear least squares:
    s = DFx\Fx;
    % [Q,R] = qr(DFx, 0);   s = R\(Q'*Fx);
    x = x-s;
    GNupdate(GNiter) = norm(s,inf);
end
lGN = log(GNupdate);
order = (lGN(3:end)-lGN(2:end-1))./(lGN(2:end-1)-lGN(1:end-2));
rate = GNupdate(2:end)./GNupdate(1:end-1);
% print norm of the update, rate and order of convergence:
GNupdate_rate_order = [GNupdate',[0;rate'],[0;0;order']]

disp('Solution: '); x,

figure('name','solution');
t_plot = linspace(t(1),t(length(t)),1000);
plot(t,m,'*', t_plot, F(x,t_plot,0),'r',...
    t_plot, F(xi,t_plot,0),'k',...
    t_plot, x(1)*exp(-x(3)*t_plot), 'r--',...
    t_plot, xi(1)*exp(-xi(3)*t_plot), 'k--','linewidth',2);
xlabel('{\bf time t}');
ylabel('{\bf Substance amounts}');
legend('data (t,m)', 'fittedPhiB','initial guess PhiB',...
     'fittedPhiA','initial guess PhiA','location','best');
print -depsc2 'EstimatedPhiAPhiB.eps';

figure('name','convergence');
semilogy(1:GNiter,GNupdate, '-o','linewidth',2);
xlabel('{\bf no. of iteration step}','fontsize',12); 
ylabel('{\bf Euclidean norm of parameter update}','fontsize',12);
print -depsc2 'GNconvergence.eps';

% CODE TEST: check if you reached your initial parameters:
% CodeTest = [xtest, x]
toc
end

function r = F(x,t,m)
% Function F to be minimized in 2-norm 
% With m=0, it corresponds to Phi_B
% Inputs:  x     4-vector,
%          t     vector/matrix of any size
%          m     vector/matrix of the size of t, or scalar
% Output:  r     vector/matrix of the size of t
r = exp(-x(4)*t) * x(2) + x(3)/(x(4)-x(3)) *...
    (exp(-x(3)*t) - exp(-x(4)*t)) * x(1) -m; 


end

function r = DF(x,t)
% Jacobian of F
% inputs:  x        4-vector,
%          t        column vector of length N
% output:  r        N x 4  matrix
d = x(4)-x(3);
ex3 = exp(-x(3)*t);
ex4 = exp(-x(4)*t);
exd = ex3-ex4;
r = [x(3)/d*exd,    ex4,...
    x(4)/d^2*x(1) * exd - x(3)*x(1)/d * t.* ex3, ...
    -x(2)*ex4 .* t + x(3)*x(1)/d *t.*ex4 - x(4)/d^2*x(1) *exd];
end
