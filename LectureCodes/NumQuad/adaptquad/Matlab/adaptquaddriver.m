function [L,X,vals,ests] = adaptquaddriver(filename,f,x,reltol,abstol)

if (nargin < 5), abstol = eps; end
if (nargin < 4), reltol = 1E-6; end
if (nargin < 3), x = 0:0.1:1; end

global res;
res.mmax = 10000;
res.level = 0;
res.X = [];
res.vals = [];
res.ests = [];

Ival = adaptquad_nc(f,x,reltol,abstol);
tp = x(1):(x(end)-x(1))/1E7:x(end);
Eval = trapz(tp,f(tp));
% Eval = quad(f,x(1),x(end),eps);

L = res.level;
X = res.X;
vals = res.vals;
ests = res.ests;

% Plotting function an points
N = 1000;
xp = x(1):(x(end)-x(1))/N:x(end);
figure('name','Adaptive quadrature');
plot3(xp,zeros(size(xp)),f(xp),'b-'); 
ax = axis;
axis([ax(1) ax(2)  0 L+1 ax(5) ax(6)]);
hold on; np = [];
for l=1:L
  plot3(X{l},l*ones(size(X{l})),zeros(size(X{l})),'r.');
  np = [np, length(X{l})];
end
grid on;
xlabel('{\bf x}','fontsize',14);
ylabel('{\bf quadrature level}','fontsize',14);
zlabel('{\bf f}','fontsize',14);
view(-220,20);

print('-depsc2',sprintf('../PICTURES/%spts.eps',filename));

figure('name','quadrature error');
semilogy(np,abs(vals-Eval),'m+',np,ests,'r*');
xlabel('{\bf no. of quadrature points}','fontsize',14);
ylabel('{\bf quadrature errors}','fontsize',14);
legend('exact error','estimated error','location','northeast');

print('-depsc2',sprintf('../PICTURES/%serrs.eps',filename));




