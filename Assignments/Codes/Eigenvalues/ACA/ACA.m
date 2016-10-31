function B = ACA(A,tol)
B = [];
[val0,i] = max(diag(A));
while true
    b = A(:,i)/sqrt(A(i,i));
    B = [B b];
    A = A - b*b';
    [val,i] = max(diag(A));
    if val < tol*val0, break; end
end
