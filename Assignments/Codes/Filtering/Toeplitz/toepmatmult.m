function y = toepmatmult(c,r,x)
    if (size(c) ~= size(r)), error('c,r length mismatch'); end
    T = toeplitz(c,r);
    y = T*x;
