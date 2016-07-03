% Setting: \Blue{$N \gg 1$}, 
% large matrix \Blue{$\VA\in\bbK^{n,n}$}
[L,U] = lu(A);
for j=1:N
    x = U\(L\b);
    b = some_function(x);
end