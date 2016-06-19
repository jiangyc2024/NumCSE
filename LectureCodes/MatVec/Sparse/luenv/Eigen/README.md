does not work with this input

	int n = 7;
	typedef Triplet<double> triplet_t;
	std::vector<triplet_t> triplets;
	
	// example in lecture document ex:env
	for(int i = 0; i < n; ++i)
		triplets.push_back( {i,i,1.0} );//diagonal
//	triplets.push_back( {0,0,1.0} );// it does not work without this line!!!
	triplets.push_back( {0,2,1.0} );
	triplets.push_back( {1,4,1.0} );
	triplets.push_back( {2,0,1.0} );
	triplets.push_back( {2,6,1.0} );
	triplets.push_back( {3,4,1.0} );
	triplets.push_back( {3,6,1.0} );
	triplets.push_back( {4,1,1.0} );
	triplets.push_back( {4,3,1.0} );
	triplets.push_back( {4,5,1.0} );
	triplets.push_back( {5,4,1.0} );
	triplets.push_back( {6,2,1.0} );
	triplets.push_back( {6,3,1.0} );
	SparseMatrix<double> A(n,n);
	A.setFromTriplets(triplets.begin(), triplets.end());
