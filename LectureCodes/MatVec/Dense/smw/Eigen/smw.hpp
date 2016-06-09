VectorXd smw(const MatrixXd &L, const MatrixXd &U, const MatrixXd &u,
				const VectorXd &v, const VectorXd &b){
	VectorXd z = U.triangularView<Upper>().solve(L.triangularView<Lower>().solve(b));
	VectorXd w = U.triangularView<Upper>().solve(L.triangularView<Lower>().solve(u));
	double alpha = 1.0 + v.dot(w);
	assert(std::abs(alpha) > std::numeric_limits<double>::epsilon()*U.lpNorm<1>()
			&& "A nearly singular");
	return z - w * v.dot(z) / alpha;
}
