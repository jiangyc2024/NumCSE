x = A.triangularView<Lower>().solve(b);				--> lower triangular
x = A.triangularView<Upper>().solve(b);				--> upper triangular
x = A.selfadjointView<Upper>().llt().solve(b);		--> Hermitian / self adjoint and positive definite
x = A.selfadjointView<Upper>().llt().solve(b);		--> Hermiatin / self adjoint (poitive or negative semidefinite)
