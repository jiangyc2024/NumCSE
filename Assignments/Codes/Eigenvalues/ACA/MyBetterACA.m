function B = MyBetterACA(A,tol)
B = [];
d = diag(A);
[val0,i] = max(d);
col = A(:,i);
while true
    b = col/sqrt(d(i));
    B = [B b];
    d = d - b.*b;
    [val,i] = max(d);
    col = A(:,i) -  B * B(i,:)';
    if val < tol*val0, break; end
end
