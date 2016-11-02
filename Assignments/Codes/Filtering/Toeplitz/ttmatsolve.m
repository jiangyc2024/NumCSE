function x = ttmatsolve(h,y)
    n = size(h,1);
    T = toeplitz(h,[h(1); zeros(n-1,1)]);
    x = T\y;
