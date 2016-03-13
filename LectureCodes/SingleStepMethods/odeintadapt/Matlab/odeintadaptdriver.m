function odeintadaptdriver(T,a,reltol,abstol)
% Simple adaptive timestepping strategy of Code~\ref{mc:odeintadapt}
% based on explicit Euler \eqref{eq:eeul} and explicit trapezoidal 
% rule \eqref{eq:exTrap}

% Default arguments
if (nargin < 4), abstol = 1E-4; end
if (nargin < 3), reltol = 1E-2; end
if (nargin < 2), a = 20; end
if (nargin < 1), T = 2; end

% autonomous ODE \Blue{$\dot{y}=\cos(ay)$} and its general solution
f = @(y) (cos(a*y).^2); sol = @(t) (atan(a*(t-1))/a);
% Initial state \Blue{$\Vy_0$}
y0 = sol(0);

% Discrete evolution operators, see Def.~\ref{def:esv}
Psilow = @(h,y) (y + h*f(y)); % Explicit Euler \eqref{eq:eeul}
% Explicit trapzoidal rule \eqref{eq:exTrap}
Psihigh = @(h,y) (y + 0.5*h*(f(y)+f(y+h*f(y)))); 

% Heuristic choice of initial timestep and \Blue{$h_{\min}$}
h0 = T/(100*(norm(f(y0))+0.1)); hmin = h0/10000;
% Main adaptive timestepping loop, see Code~\ref{mc:odeintadapt}
[t,y,rej,ee] = odeintadapt_ext(Psilow,Psihigh,T,y0,h0,reltol,abstol,hmin);

% Plotting the exact the approximate solutions and rejected timesteps
figure('name','solutions');
tp = 0:T/1000:T; plot(tp,sol(tp),'g-','linewidth',2); hold on;
plot(t,y,'r.'); 
plot(rej,0,'m+');
title(sprintf('Adaptive timestepping, rtol = %f, atol = %f, a = %f',reltol,abstol,a));
xlabel('{\bf t}','fontsize',14);
ylabel('{\bf y}','fontsize',14);
legend('y(t)','y_k','rejection','location','northwest');
print -depsc2 '../PICTURES/odeintadaptsol.eps';

fprintf('%d timesteps, %d rejected timesteps\n',length(t)-1,length(rej));

% Plotting estimated and true errors
figure('name','(estimated) error');
plot(t,abs(sol(t) - y),'r+',t,ee,'m*');
xlabel('{\bf t}','fontsize',14);
ylabel('{\bf error}','fontsize',14);
legend('true error |y(t_k)-y_k|','estimated error EST_k','location','northwest');
title(sprintf('Adaptive timestepping, rtol = %f, atol = %f, a = %f',reltol,abstol,a));
print -depsc2 '../PICTURES/odeintadapterr.eps';



