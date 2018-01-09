function gamma = data_fit_normal(t,f)

t = t(:);
f = f(:);

A = make_A(t);

AtA = A'*A;
Atf = A'*f;

gamma = AtA\Atf;