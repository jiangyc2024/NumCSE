function x = arrowsys_slow(d,c,b,alpha,y)
A = [diag(d),c;transpose(b),alpha];
x = A\y;
