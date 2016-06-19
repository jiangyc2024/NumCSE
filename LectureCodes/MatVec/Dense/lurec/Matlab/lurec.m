function A = lurec(A)
% insitu recursive LU-factorization
if(size(A,1)>1)
 fac = A(2:end,1)/A(1,1);
 C = lurec(A(2:end,2:end)...
 -fac*A(1,2:end));
 A=[A(1,:);fac,C];
end
