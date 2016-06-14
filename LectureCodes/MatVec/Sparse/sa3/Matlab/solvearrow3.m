function x = solvearrow3(alpha,b,c,d,y)
A = [alpha, b'; c, spdiags(d,0,length(d),length(d))];
x = A\y;