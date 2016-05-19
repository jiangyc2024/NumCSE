function y = lrtrimult(A,B,x)
y = triu(A*B')*x;
