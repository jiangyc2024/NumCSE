function y = toepmult(c,r,x)
    n = size(c,1);
    if (n ~= size(r,1)), error('c,r length mismatch'); end
    y = pconvfft([c;0;r(n:-1:2)], [x;zeros(n,1)]);
    y = y(1:n);
