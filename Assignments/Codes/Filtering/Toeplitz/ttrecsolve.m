function x = ttrecsolve(h,y,l)
    if (l == 0) 
        x(1) = y(1)/h(1);
    else
        n = 2^l; m = n/2;
        if ((size(h,1) ~= n) || (size(y,1) ~= n))
            error('Size mismatch'); end
        x1 = ttrecsolve(h(1:m), y(1:m),l-1);
        y2 = y(m+1:n) - toepmult(h(m+1:n), h(m+1:-1:2),x1);
        x2 = ttrecsolve(h(1:m,y2,l-1);
        x  = [x1; x2];
    end
