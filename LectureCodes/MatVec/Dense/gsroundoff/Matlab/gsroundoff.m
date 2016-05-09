% MATLAB script demonstrating the effect of roundoff on the result of Gram-Schmidt orthogonalization
format short; % Print only a few digits in outputs
% Create special matrix the so-called Hilbert matrix: $\MAc{i}{j} = (i+j-1)^{-1}$
A = hilb(10);       % 10x10 Hilbert matrix 
Q = gramschmidt(A); % Gram-Schmidt orthogonalization of columns of A
% \Magenta{Test orthonormality} of column of Q, which should be an \cor{orthogonal}
% matrix according to theory
I = Q'*Q,  % Should be the unit matrix, but isn't !
 
% MATLAB's internal Gram-Schmidt orthogonalization by \cor{QR-decomposition}
[Q1,R1] = qr(A,0); 
D = A - Q1*R1,  % Check whether we get the expected result
I1 = Q1'*Q1,    % Test orthonormality