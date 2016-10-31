% initialization for rank 1 circuit modification
clear all;
Res   = 50;
C     =  3;
L     =  1;
w     = 0.2;              % angular frequency
A     = r1circuit_getA(Res,Res,C,L,w);
[Q,R] = qr(full(A));      % full in order to use qrupdate

U_in  = 1;                % input voltage
b     = sparse(16,1);     % right hand side
b(7)  = U_in / Res;

U_all   = R\(Q'*b);
I_15_16 = (U_all(15)-U_all(16))/Res