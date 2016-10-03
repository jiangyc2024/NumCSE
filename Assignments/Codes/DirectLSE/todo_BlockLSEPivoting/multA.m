function x = multA (d1,c,d2,x)
n = length(x)/2;
x = [d1; d2] .* x + [c; c] .* [x(n+1:2*n); x(1:n)];
