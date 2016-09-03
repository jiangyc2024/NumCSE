#include <Eigen/Dense>
#include <Eigen/Sparse>

template <class SimilarityFunction>
std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> 
	imgsegmat(Eigen::MatrixXd P, SimilarityFunction Sfun)
{
	using Triplet = Eigen::Triplet<double>;

	int n = P.rows();
	int m = P.cols();
	int N = n*m;

	std::vector<Triplet> WTriplets;
	std::vector<Triplet> DTriplets;
	WTriplets.reserve(4*N);
	DTriplets.reserve(4*N);

	int k = 1;
	for (int ni=0; ni<n; ++ni)
	{
		for (int mi=0; mi<m; ++mi)
		{
			int mni = mi*n+ni; //global index

			// call Sfun with adjacent pixels
			if (ni-1>=0) // bottom pixel
			{
				double value = Sfun(P(ni,mi),P(ni-1,mi));
				WTriplets.push_back(Triplet(mni, mni-1, value));
				DTriplets.push_back(Triplet(mni, mni, value));
				k = k + 1;
			}
			
			if (ni+1<n) // top pixel
			{
				double value = Sfun(P(ni,mi),P(ni+1,mi));
				WTriplets.push_back(Triplet(mni, mni+1, value));
				DTriplets.push_back(Triplet(mni, mni, value));
				k = k + 1;
			}
			
			if (mi-1>=0) // left pixel
			{
				double value = Sfun(P(ni,mi),P(ni,mi-1));
				WTriplets.push_back(Triplet(mni, mni-n, value));
				DTriplets.push_back(Triplet(mni, mni, value));
				k = k + 1;
			}

			if (mi+1<m) // right pixel
			{
				double value = Sfun(P(ni,mi),P(ni,mi+1));
				WTriplets.push_back(Triplet(mni, mni+n, value));
				DTriplets.push_back(Triplet(mni, mni, value));
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

