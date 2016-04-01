function z = zerosquadpolstab(alpha,beta)
% MATLAB function computing the real zeros of a quadratic polynomial
% $\xi\to \xi^2+\alpha\xi+\beta$ by means of the familiar discriminant
% formula $\xi_{1,2} = \frac{1}{2}(-\alpha\pm\sqrt{\alpha^2-4\beta})$ and
% This is a \emph{stable} implementation based on \emph{Vieta's theorem}. 
D = alpha^2-4*beta; % discriminant
if (D < 0), z = []; 
else
    wD = sqrt(D); 
    % Use discriminant formula only for zero far away from $0$ 
    % in order to \com{avoid cancellation}. For the other zero
    % use Vieta's formula. 
    if (alpha >= 0)
        t = 0.5*(-alpha-wD); % \Label[line]{zqs:1}
        z = [t;beta/t];
    else
        t = 0.5*(-alpha+wD); % \Label[line]{zqs:2}
        z = [beta/t;t];
    end
end

