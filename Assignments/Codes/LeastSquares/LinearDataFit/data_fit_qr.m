function gamma = data_fit_qr(t,f)

t = t(:);
f = f(:);

A = make_A(t);

[Q,R] = qr(A,0);
gamma = R\(Q'*f);