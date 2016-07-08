function x1 = ModNewtStep (x0, F, DF)
y0 = x0 + F(x0) / DF(x0);    % version for scalar eqs only
x1 = y0 - F(y0) / DF(x0);

% y0 = x0 + DF(x0) \ F(x0);  % version for scalar and vector eqs
% x1 = y0 - DF(x0) \ F(y0);
