function x = ModNewtSys( A, c, tol )
c = c(:);         % make sure it is column
F = @(x) A*x + c .* exp(x);
DF = @(x) A + diag( c .* exp(x) );
x0 = zeros( size(c) );

x = x0;
it = 0;
RelIncr = 1;
res = norm(F(x));
while ( RelIncr > tol  &&  it < 100 )
    y = x(:,end) + DF(x(:,end)) \ F(x(:,end));
    x_new = y - DF(x(:,end)) \ F(y);
    RelIncr = norm( x(:,end) - x_new ) / norm(x_new);
    x = [ x, x_new ];
    it = it + 1;
    res = [res, norm(F(x(:,end))) ];
end
disp('Sequence of residuals:')
res

% return only final values and forget history
x = x(:, end);

% test:
%close all; plot(ModNewtSys(gallery('poisson',20),(1:400),1e-10));
