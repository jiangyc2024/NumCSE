function A = blockgs(A)
%in-situ recursive Gaussian elimination, no pivoting
%right hand side in rightmost column of \Blue{$\VA$}: A(:,end)
n=size(A,1); 
if(n>1)
 C=blockgs(A(2:end,2:end)-A(2:end,1)...
     *A(1,2:end)/A(1,1));
 A=[A(1,:);zeros(n-1,1),C];
end
