% Example: numerical instability without pivoting
A = [5.0E-17 , 1; 1 , 1];
b = [1;2];
x1 = A\b,
x2 =gausselim(A,b), % see Code~\ref{gausselim}
[L,U] = lufak(A); % see Code~\ref{mc:lufak}
z = L\b; x3 = U\z,
