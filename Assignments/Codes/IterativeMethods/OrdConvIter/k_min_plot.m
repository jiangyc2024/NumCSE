% MATLAB Plotting script for Course "Mumerical Methods for CS(E)", Order-p convergent
% iteration

C               = 2; % Costant 
p               = 1.5; % Order p 
eps_max         = C^(1/(1-p)); % Threshold for intial error 

ngp             = 100; % number of grid points
eps_lin         = linspace(0, eps_max, ngp);
tau_lin         = linspace(-6, -1, ngp); % Logarithm of error reduction 
[eps_msh,tau_msh] = meshgrid(eps_lin(2:(end-1)),tau_lin(2:(end-1)));

kmin = @(eps, C, p, tau) ceil(log( ( log(tau) + (1/(p-1)).*log(C) ) ./ log(C^(1/(p-1)) .* eps)  ) ./ log(p) );
k = kmin(eps_msh, C, p, 10.^tau_msh);

% Consider only gridpoints where: eps larger as tau
for ne = 1:ngp-2
    for nt = 1:ngp-2
        if (eps_lin(ne+1) < 10^tau_lin(nt+1))
           k(nt,ne) = 0;
        end
    end
end

% Plotting
figure('Name','k_min plot');
mesh(eps_msh,tau_msh,k);
title('Minimal number of iterations for error < \tau');
xlabel('\epsilon_0');
ylabel('log_{10}\tau');
zlabel('k_{min}');
xlim([0,eps_max]);
ylim([-6,-1]);

% % Just a color plot
% figure('Name','k_min plot');
% pcolor(eps_msh,tau_msh,k)
% colorbar()
% title('Minimal number of iterations for error < \tau')
% xlabel('\epsilon_0')
% ylabel('\tau')
% xlim([0,eps_max])
% ylim([-6,-1])
% shading flat