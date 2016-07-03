



std::pair<Eigen::MatrixXd, Eigen::MatrixXd> imgsetmat(Eigen::MatrixXd P)
{
	int n = P.rows();
	int m = P.cols();

	std::vector<Eigen::Triplet<double>> triplets(4*n*m);
	spdata = zeros(4*n*m,1); spidx  = zeros(4*n*m,2);
	int k = 1;
	for (int ni=0; ni<=n; ++ni)
	{
		for (int mi=0; mi<=m; ++mi)
		{
			int mni = mi*n+ni;
			if (ni-1>=0) 
			{
				double value = Sfun(P(ni,mi),P(ni-1,mi));
				triplets.push_back(Eigen::Triplet(mni, mni-1));
				k = k + 1;
			}
			
			if (ni+1<n)
			{
				double value = Sfun(P(ni,mi),P(ni+1,mi));
				triplets.push_back(Eigen::Triplet(mni, mni+1));
				k = k + 1;
			}
			
			if (mi-1>0)
			{
				spidx(k,:) = [mni,mni-n];
				spdata(k) = Sfun(P(ni,mi),P(ni,mi-1));
				k = k + 1;
			}

			if (mi+1<=m), spidx(k,:) = [mni,mni+n];
			spdata(k) = Sfun(P(ni,mi),P(ni,mi+1));
			k = k + 1;
			end
		}
	}

	// Efficient initialization, see Sect.~\ref{sec:spml}, Ex.~\ref{ex:spinit}
	Eigen::SparseMatrix<double> W(n*m,n*m);
	W.setFromTriplets(triplets);
	W.compress();

	W = sparse(spidx(1:k-1,1),spidx(1:k-1,2),spdata(1:k-1),n*m,n*m);
	D = spdiags(full(sum(W')'),[0],n*m,n*m);
	A = D-W;
 
}

