#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace imgsegmat {


template <class SimilarityFunction>
inline std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> 
	imgsegmat(Eigen::MatrixXd P, SimilarityFunction Sfun)
{
	using Triplet = Eigen::Triplet<double>;
	using index_t = Eigen::Index;

	const index_t n = P.rows();
	const index_t m = P.cols();
	const index_t N = n*m;

	std::vector<Triplet> WTriplets;
	std::vector<Triplet> DTriplets;
	WTriplets.reserve(4*N);
	DTriplets.reserve(4*N);

	index_t k = 1;
	for (index_t ni=0; ni<n; ++ni)
	{
		for (index_t mi=0; mi<m; ++mi)
		{
			const index_t mni = mi*n+ni; //global index

			// call Sfun with adjacent pixels
			if (ni-1>=0) // bottom pixel
			{
				const double value = Sfun(P(ni,mi),P(ni-1,mi));
				WTriplets.emplace_back(mni, mni-1, value);
				DTriplets.emplace_back(mni, mni, value);
				k = k + 1;
			}
			
			if (ni+1<n) // top pixel
			{
				const double value = Sfun(P(ni,mi),P(ni+1,mi));
				WTriplets.emplace_back(mni, mni+1, value);
				DTriplets.emplace_back(mni, mni, value);
				k = k + 1;
			}
			
			if (mi-1>=0) // left pixel
			{
				const double value = Sfun(P(ni,mi),P(ni,mi-1));
				WTriplets.emplace_back(mni, mni-n, value);
				DTriplets.emplace_back(mni, mni, value);
				k = k + 1;
			}

			if (mi+1<m) // right pixel
			{
				const double value = Sfun(P(ni,mi),P(ni,mi+1));
				WTriplets.emplace_back(mni, mni+n, value);
				DTriplets.emplace_back(mni, mni, value);
				k = k + 1;
			}
		}
	}


	// Efficient initialization, see Sect.~\ref{sec:spml}, Ex.~\ref{ex:spinit}
	Eigen::SparseMatrix<double> W(N,N);
	Eigen::SparseMatrix<double> D(N,N);
	Eigen::SparseMatrix<double> A(N,N);

	W.setFromTriplets(WTriplets.begin(), WTriplets.end());
	D.setFromTriplets(DTriplets.begin(), DTriplets.end());
	W.makeCompressed();
	D.makeCompressed();

	A = D-W;
	return std::make_pair(A,D);
}


} //namespace imgsegmat