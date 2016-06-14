% Example: importance of scale-invariant pivoting
epsilon = 5.0E-17;
A = [epsilon , 1; 1 , 1]; b = [1;2];
D = [1/epsilon, 0; 0 ,1];
A = D*A; b = D*b;
x1 = A\b, % MATLAB internal Gaussian elimination
x2 =gausselim(A,b), % see Code~\ref{gausselim}
[L,U] = lufak(A); % see Code~\ref{mc:lufak}
z = L\b; x3 = U\z,
