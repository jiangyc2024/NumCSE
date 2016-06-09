% Setting: \Blue{$N \gg 1$}, 
% large matrix \Blue{$\VA\in\bbK^{n,n}$}
for j=1:N
    x = A\b;
    b = some_function(x);
end