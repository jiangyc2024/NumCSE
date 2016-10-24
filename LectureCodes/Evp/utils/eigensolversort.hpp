#include <Eigen/Dense>
#include <vector>


// sorts eigenvalues and eigenvectors from eigensolver ev
std::pair<Eigen::VectorXd, Eigen::MatrixXd> eigensolversort(Eigen::EigenSolver<Eigen::MatrixXd> ev, bool sortAbs = false)
{
	Eigen::MatrixXd eigenvectors = ev.eigenvectors().real();
	Eigen::VectorXd eigenvalues = ev.eigenvalues().real();
	if (sortAbs) eigenvalues = ev.eigenvalues().cwiseAbs();

	int n = eigenvalues.size();
	std::vector<std::pair<double, Eigen::VectorXd>> res;
	res.reserve(n);

	for (int i=0; i<n; ++i)
	{
		res.push_back(std::make_pair(eigenvalues(i), Eigen::VectorXd(eigenvectors.col(i))));
	}

	std::sort(res.begin(), res.end(),
		[](const std::pair<double, Eigen::VectorXd> &a, const std::pair<double, Eigen::VectorXd> &b)
		{ return a.first < b.first; });

	Eigen::VectorXd eigenval(n);
	Eigen::MatrixXd eigenvecs(n,n);
	
	for (int i=0; i<n; ++i)
	{
		eigenval(i) = res[i].first;
		eigenvecs.col(i) = res[i].second;
	}

	return std::make_pair(eigenval, eigenvecs);
}
