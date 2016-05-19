% MATLAB script for testing the computation of the zeros of a parabola
res = []; % array for storing results of computations. 
gammavec = 2:10:1000; % Define test cases
for gamma=gammavec
    % Polynomial $p(\xi) = (\xi-\gamma)(\xi-\frac{1}{\gamma})$
   alpha = -(gamma + 1/gamma); % compute coefficients of polynomial
   beta = 1.0; 
   z1 = zerosquadpol(alpha,beta); % \emph{Unstable} way to compute zeros
   z2 = zerosquadpolstab(alpha,beta); % \emph{Stable} implementation
   % Compute relative errors of the left zero
   ztrue = 1/gamma; z2true = gamma;
   res = [res; gamma, abs((z1(1)-ztrue)/ztrue),...
       abs((z2(1)-ztrue)/ztrue),abs((z1(2)-z2true)/z2true)];    
end
% Graphical output of relative error of roots computed by unstable
% implementation
figure; plot(res(:,1),res(:,2),'r+',res(:,1),res(:,4),'b*');  
legend('small root','large root','location','best');
xlabel('{\bf \gamma}','fontsize',14);
ylabel('{\bf relative errors in \xi_1, \xi_2}','fontsize',14);
title('Roots of a parabola computed in an unstable manner');   
print -depsc2 '../PICTURES/zqperrinstab.eps';
% Graphical output of relative errors (comparison)
figure; plot(res(:,1),res(:,2),'r+',res(:,1),res(:,3),'m*');  
legend('unstable','stable');
xlabel('{\bf \gamma}','fontsize',14);
ylabel('{\bf relative error in \xi_1}','fontsize',14);
title('Roundoff in the computation of zeros of a parabola');
print -depsc2 '../PICTURES/zqperr.eps';

