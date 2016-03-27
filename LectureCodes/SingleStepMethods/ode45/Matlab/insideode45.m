function varargout = ode45(ode,tspan,y0,options,varargin)
% Processing of input parameters omitted
% $\vdots$
% Initialize method parameters, \emph{c.f.} \Red{Butcher scheme} \eqref{eq:BSexpl}
pow = 1/5;
A = [1/5, 3/10, 4/5, 8/9, 1, 1];
B = [
    1/5         3/40    44/45   19372/6561      9017/3168       35/384
    0           9/40    -56/15  -25360/2187     -355/33         0
    0           0       32/9    64448/6561      46732/5247      500/1113
    0           0       0       -212/729        49/176          125/192
    0           0       0       0               -5103/18656     -2187/6784
    0           0       0       0               0               11/84
    0           0       0       0               0               0
    ];
E = [71/57600; 0; -71/16695; 71/1920; -17253/339200; 22/525; -1/40];
% $\vdots$  (choice of stepsize and main loop omitted)
% ADVANCING ONE STEP.
hA = h * A;
hB = h * B;
f(:,2) = feval(odeFcn,t+hA(1),y+f*hB(:,1),odeArgs{:});
f(:,3) = feval(odeFcn,t+hA(2),y+f*hB(:,2),odeArgs{:});
f(:,4) = feval(odeFcn,t+hA(3),y+f*hB(:,3),odeArgs{:});
f(:,5) = feval(odeFcn,t+hA(4),y+f*hB(:,4),odeArgs{:});
f(:,6) = feval(odeFcn,t+hA(5),y+f*hB(:,5),odeArgs{:});

tnew = t + hA(6);
if done, tnew = tfinal; end  % Hit end point exactly.
h = tnew - t;      % Purify h.    
ynew = y + f*hB(:,6);  
% $\vdots$ (stepsize control, see Sect.~\ref{sec:ssctrl} dropped  