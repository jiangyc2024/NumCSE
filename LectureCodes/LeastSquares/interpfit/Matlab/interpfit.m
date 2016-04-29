% Comparison of polynomial interpolation and polynomial fitting
% (``Quick and dirty'' MATLAB implementation, see \ref{subsec:interp-algorithms})
f = @(x) 1./(1+x.^2); % Function providing data points
d = 10;              % Polynomial degree
tip = -5+(0:d)*10/d;       % \Blue{$d+1$} nodes for interpolation 
tft = -5+(0:3*d)*10/(3*d); % \Blue{$3d+1$} Nodes for polynomial fitting 
pip = polyfit(tip,f(tip),d); % Interpolating polynomial (degree \Blue{$d$})
pft = polyfit(tft,f(tft),d); % Fitting polynomial (degree \Blue{$d$})
x = -5+(0:1000)/100;
h = plot(x,f(x),'g--',...
         x,polyval(pip,x),'b-',...
         x,polyval(pft,x),'r-',...
         tip,f(tip),'b*');
set(h(1),'linewidth',2);
xlabel('{\bf t}','fontsize',14);
ylabel('{\bf y}','fontsize',14);
legend('function f','interpolating polynomial','fitting polynomial','location','north');
print -depsc2 '../PICTURES/interpfit.eps';


