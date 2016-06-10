template <class solver_t>
VectorXd arrowsys_sparse(const VectorXd &d, const VectorXd &c, 
		const VectorXd &b, const double alpha, const VectorXd &y){
	int n = d.size();
	SparseMatrix<double> A(n+1, n+1); // default: column-major
	VectorXi reserveVec = VectorXi::Constant(n+1, 2); // nnz per col
	reserveVec(n) = n+1;		// last full col
	A.reserve(reserveVec);
	for(int j = 0; j < n; ++j){	// initalize along cols for efficency
		A.insert(j,j) = d(j);	// diagonal entries
		A.insert(n,j) = b(j);	// bottom row entries
	}
	for(int i = 0; i < n; ++i){
		A.insert(i,n) = c(i);	// last col
	}
	A.insert(n,n) = alpha;		// bottomRight entry
	A.makeCompressed();
	return solver_t(A).solve(y);
}
