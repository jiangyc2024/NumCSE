function z = zerosquadpol(alpha,beta)
% MATLAB function computing the zeros of a quadratic polynomial
% $\xi\to \xi^2+\alpha\xi+\beta$ by means of the familiar discriminant
% formula $\xi_{1,2} = \frac{1}{2}(-\alpha\pm\sqrt{\alpha^2-4\beta})$. However
% this implementation is \emph{vulnerable to round-off}! The zeros are
% returned in a column vector
D = alpha^2-4*beta; % discriminant
if (D < 0), z = []; % No real zeros
else
    % The famous discriminant formula
    wD = sqrt(D); 
    z = 0.5*[-alpha-wD;-alpha+wD];
end
