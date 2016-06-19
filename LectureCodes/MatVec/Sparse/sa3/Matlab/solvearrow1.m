function x = solvearrow1(alpha,b,c,d,y)
A = [alpha, b'; c, diag(d)];
x = A\y;