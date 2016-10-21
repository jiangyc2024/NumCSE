% INPUT
% t, y - parameters, such that f(t(i)) = y(i)
% OUTPUT
% alpha, beta - values such that alpha*exp(beta*t) approximates f, through
% least squares.

function [alpha, beta] = ExpoFuncFit(t, y)
    % Compute the solution of the linearised problem
    x = linearregression(t, log(y));
    % Set the proper output values
    alpha = exp(x(2));
    beta = x(1);
end

function x = linearregression(t,y)
% Solution of linear regression problem (fitting of a line to data) for 
% data points (t_i,y_i), i=1,\ldots,n, passed in the 
% column vectors t and y. 
% The return value is a 2-vector, containing the slope of the fitted line 
% in x(1) and its offset in x(2)
n = length(t); 
if (length(y) ~= n)
    error('data size mismatch'); 
end
% Coefficient matrix of \textbf{overdetermined linear system}
A = [t,ones(n,1)]; 
% Determine least squares solution by using MATLAB's ``backslash''
% operator}
x = A\y;
end
