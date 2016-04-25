function chemstiff
% Simulation of kinetics of coupled chemical reactions with vastly different reaction
% rates, see \eqref{eq:chemstiff} for the ODE model.
% reaction rates \Blue{$k_1,k_2,k_3,k_4$}, \Magenta{$k_1,k_2 \gg k_3,k_4$}.
k1 = 1E4; k2 = 1E3; k3 = 10; k4 = 1;
% definition of right hand side function for ODE solver
fun = @(t,y) ([
		-k1*y(1)*y(2) + k2*y(3) - k3*y(1)*y(3) +  k4*y(4);
		-k1*y(1)*y(2) + k2*y(3);
		k1*y(1)*y(2) - k2*y(3) - k3*y(1)*y(3) + k4*y(4);
		k3*y(1)*y(3) - k4*y(4)]);

tspan = [0 1]; % Integration time interval
L = tspan(2)-tspan(1); % Duration of simulation
y0 = [1;1;10;0]; % Initial value \Blue{$\Vy_0$}
% compute ``exact'' solution, using \texttt{ode113} with tight error tolerances
options = odeset('reltol',10*eps,'abstol',eps,'stats','on');
% get the 'exact' solution using ode113
[tex,yex] = ode113(fun,[0 1],y0,options);
