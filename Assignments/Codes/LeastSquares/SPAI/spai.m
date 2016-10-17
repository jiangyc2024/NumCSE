function B = spai(A)
B = A;
for i=1:size(A,1)
    idx = find(A(:,i));
    B(idx,i) = ( A(:,idx)' * A(:,idx) ) \ A(i,idx)';
end
