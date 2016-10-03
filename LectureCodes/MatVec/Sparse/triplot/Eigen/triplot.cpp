#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>

#include <figure/figure.hpp>

#include <mgl2/mgl.h>

void TriPlot(Eigen::Matrix<double, -1, -1, Eigen::RowMajor> &T, const Eigen::VectorXd &x, const Eigen::VectorXd &y){
	mglGraph gr;
	mglData xd(x.data(), x.size()),
			yd(y.data(), y.size()),
			Td(T.rows(), T.cols(), T.data());
	gr.SetRanges(0,1,0,1);
	//gr.Grid("","h="); 
	gr.TriPlot(Td, xd, yd,"#b");
	gr.Plot(xd, yd, " r*");
	gr.Label(xd, yd,"%n"); // give each point an individual number
	gr.Axis();
	gr.WriteEPS("meshplot_cpp.eps");
	
}
void TriPlot(Eigen::MatrixXd &T, const Eigen::VectorXd &x, const Eigen::VectorXd &y){
	// NOTE: MathGL's mglData constructor needs the matrix in RowMajor!
	Eigen::Matrix<double, -1, -1, Eigen::RowMajor> TRow(T);
	TriPlot(TRow,x,y);
}


#include <vector>
void processmesh(const Eigen::MatrixXi &T, Eigen::MatrixXi & E, Eigen::MatrixXi &Eb){
	
	// Number of nodes of the triangular mesh
	int N = T.maxCoeff()+1;
	// Number of triangles of the mesh 
	int M = T.rows();
	
	//std::cout << "T\n" << T << std::endl;
	
	std::vector<Eigen::Triplet< int> > triplets;
	triplets.reserve(3*M);
	// loop over all triangles
	for(int k = 0; k < M; ++k){
		// loop over all possible  edge combinations
		for(int i = 0; i < 2; ++i){
			for(int j = i+1; j < 3; ++j){
				if(T(k,i) < T(k,j))	// insert combinations sorted
					triplets.push_back({T(k,i) ,T(k,j), 1});
				else
					triplets.push_back({T(k,j) ,T(k,i), 1});
			}
		}
	}
	
	// build sparse Matrix
	Eigen::SparseMatrix<int> A(N,N);
	//std::cout << N << std::endl;
	/*for(auto it:triplets){
		std::cout << "{" << it.row() <<"," << it.col() << "," << it.value() << "}\n";
	}*/
	A.setFromTriplets(triplets.begin(), triplets.end());
	//std::cout << A << std::endl;
	
	// Iteration over non-zero elements
	std::vector<int> e;
	std::vector<int> eb;
	// works only for col major
	for(int i = 0; i < A.outerSize(); ++i){
		for(Eigen::SparseMatrix<int>::InnerIterator it(A, i); it; ++it){
			if(it.value() == 1){	// boundary edge
				eb.push_back(it.row());
				eb.push_back(it.col());
			}
			e.push_back(it.row());
			e.push_back(it.col());
		}
	}
	E = Eigen::Matrix<int, -1, -1,Eigen::RowMajor>::Map(e.data(), e.size()/2, 2);
	Eb = Eigen::Matrix<int, -1, -1,Eigen::RowMajor>::Map(eb.data(), eb.size()/2, 2); // still transformation to rowmajor
	//std::cout << "E\n" << E << std::endl<< std::endl << "Eb\n" << Eb << std::endl;
	
}

void getinfo(const Eigen::MatrixXi &T, const Eigen::MatrixXi & E, Eigen::MatrixXi &ET){
	// Number of edges
	int L = E.rows();
	// Number of nodes of the triangular mesh
	int N = T.maxCoeff()+1;
	
	std::vector<Eigen::Triplet< int> > triplets;
	triplets.reserve(L);
	for(int i = 0; i < L; ++i){
		triplets.push_back({E(i,0), E(i,1), i});
		triplets.push_back({E(i,1), E(i,0), i});	// symmetrical
	}
	// build sparse Matrix
	Eigen::SparseMatrix<int> A(N,N);
	A.setFromTriplets(triplets.begin(), triplets.end());
	
	std::vector<int> ET_data;
	for(int i = 0; i < T.rows(); ++i){
		ET_data.push_back(A.coeff( T(i,1),T(i,2) ));
		ET_data.push_back(A.coeff( T(i,0),T(i,2) ));
		ET_data.push_back(A.coeff( T(i,0),T(i,1) ));
		
		
	}
	ET = Eigen::Matrix<int, -1, -1,Eigen::RowMajor>::Map(ET_data.data(), ET_data.size()/3, 3); // still transformation to rowmajor
	//std::cout << "ET\n" << ET << std::endl;
	//std::cout << "size" << N << std::endl;
}


void refinemesh(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const Eigen::MatrixXi& T, Eigen::VectorXd& x_ref, Eigen::VectorXd& y_ref, Eigen::MatrixXi &T_ref){
	Eigen::MatrixXi E, Eb, ET;
	processmesh(T,E,Eb);
	
	x_ref.resize(x.size() + E.rows());
	y_ref.resize(y.size() + E.rows());
	x_ref.head(x.size()) = x;
	y_ref.head(y.size()) = y;
	for(int i = 0; i < E.rows(); ++i){
		x_ref(i+x.size()) = ( x(E(i,0)) + x(E(i,1)) )/2.0;
		y_ref(i+y.size()) = ( y(E(i,0)) + y(E(i,1)) )/2.0;
	}
	
	getinfo(T, E, ET);
	
	assert(T.rows() == ET.rows());
	
	// Build a new list of triangles
	int Nt = T.rows();
	int Nv = x.size();
	T_ref.resize(4*Nt, 3);
	/*
	std::cout << "T=" << std::endl;
	std::cout << T << std::endl;
	std::cout << "ET=" << std::endl;
	std::cout << ET << std::endl;
	*/
	for(int i = 0; i < Nt; ++i){
		T_ref(4*i,0) = T(i,0); T_ref(4*i,1) = ET(i,1) + Nv; T_ref(4*i,2) = ET(i,2) + Nv;
		T_ref(4*i+1,0) = T(i,1); T_ref(4*i+1,1) = ET(i,0) + Nv; T_ref(4*i+1,2) = ET(i,2) + Nv;
		T_ref(4*i+2,0) = T(i,2); T_ref(4*i+2,1) = ET(i,0) + Nv; T_ref(4*i+2,2) = ET(i,1) + Nv;
		T_ref(4*i+3,0) = ET(i,0) + Nv; T_ref(4*i+3,1) = ET(i,1) + Nv; T_ref(4*i+3,2) = ET(i,2) + Nv;
		
	}
	std::cout << "T_ref=\n" << T_ref << std::endl;
}


#include <algorithm>
#include <iterator>
void smoothmesh(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const Eigen::MatrixXi& T, Eigen::VectorXd& xs, Eigen::VectorXd& ys){

	// Number of nodes of the mesh
	int Nv = x.size();
	Eigen::MatrixXi E, Eb;
	processmesh(T,E,Eb);								//<<<<--------------- processmesh hat ein fehler
	int Ne = E.rows();
	std::cout << "Eb\n" << Eb << std::endl;
	std::vector<int> bd_nodes(Eb.data(), Eb.data() + Eb.size());
	std::sort(bd_nodes.begin(), bd_nodes.end());
	auto last = std::unique(bd_nodes.begin(), bd_nodes.end());
	bd_nodes.erase(last, bd_nodes.end());
	int Nb = bd_nodes.size();
	
	std::cout << "bd Nodes: ";
	for(auto bdIdx:bd_nodes){
		std::cout << bdIdx << ", ";
	}
	std::cout << std::endl;
	
	// get boundary nodes
	Eigen::VectorXd x_bd(Nb ), y_bd(Nb);
	int counter = 0;
	for(auto idx : bd_nodes){
		x_bd(counter) = x(idx);
		y_bd(counter) = y(idx);
		++counter;
	}
	
	
	// get interior nodes
	std::vector<int> int_nodes;
	int Ni = Nv-Nb;
	int_nodes.reserve(Ni);
	
	std::vector<int> nodes_map_inv(Nv);
	std::iota(nodes_map_inv.begin(), nodes_map_inv.end(), 0);
	
	std::set_difference(nodes_map_inv.begin(), nodes_map_inv.end(), bd_nodes.begin(), bd_nodes.end(), std::inserter(int_nodes, int_nodes.begin()));
	/*
	std::cout << "int Nodes: ";
	for(auto intIdx:int_nodes){
		std::cout << intIdx << ", ";
	}
	std::cout << std::endl;
	*/
	// map from old positions to new position in matrix
	// [i] returns new position of i(old position)
	std::vector<int> nodes_map(int_nodes.begin(), int_nodes.end());
	nodes_map.insert(nodes_map.end(), bd_nodes.begin(), bd_nodes.end());
	
	// inverse map to original node position

	std::sort(nodes_map_inv.begin(), nodes_map_inv.end(),
				[&nodes_map](size_t a, size_t b){ return nodes_map[a] < nodes_map[b];});
	
	
	/*
	Eigen::VectorXd x_int(Nv - bd_nodes.size()), y_int(Nv - bd_nodes.size());
	counter = 0;
	for(auto idx : int_nodes){
		x_int(counter) = x(idx);
		y_int(counter) = y(idx);
		++counter;
	}
	*/
	
	
	
	// Counting neighbours and insert $-1$ entries
	std::vector<Eigen::Triplet<int> > triplets_int;
	std::vector<Eigen::Triplet<int> > triplets_bd;
	// reserve ??
	Eigen::VectorXi neighbours(Nv);
	neighbours.setZero();
	//std::cout << "E\n" << E << std::endl;
	
	for(int i = 0; i < E.rows(); ++i){
		
		int idx1 = nodes_map_inv[E(i,0)];
		int idx2 = nodes_map_inv[E(i,1)];
		if(idx1 < Ni && idx2 < Ni){	// edge belongs to interior area
			triplets_int.push_back({idx1, idx2, -1});
			triplets_int.push_back({idx2, idx1, -1});
		}
		else if(idx1 < Ni){	// only 1st edge node belongs to interior area
			triplets_bd.push_back({idx1, idx2-Ni, -1});
		}
		else if(idx2 < Ni){	// only 2nd edge node belongs to interior area
			triplets_bd.push_back({idx2, idx1-Ni, -1});

		}
		++neighbours(idx1);
		++neighbours(idx2);
	}
	
	//std::cout << "neighbours\n" << neighbours << std::endl;
	
	for(int i = 0; i < Ni; ++i){	// interior
		triplets_int.push_back({i,i,neighbours(i)});
	}
	for(int i = Ni; i < std::min(Ni + Nb, Ni + Ni); ++i){	// boundary
		//std::cout << i-Ni << std::endl;
		//triplets_bd.push_back({i-Ni,i-Ni,neighbours(i)});		Warum ist das so??? (auskommentiert ist richtig)
	}

	//std::cout << "x_bd\n" << x_bd << std::endl;
	//std::cout << "y_bd\n" << y_bd << std::endl;


	Eigen::SparseMatrix<double> A_int(Ni,Ni);
	A_int.setFromTriplets(triplets_int.begin(), triplets_int.end());
	//std::cout << "test " << Ni << " " << Nb << std::endl;
	Eigen::SparseMatrix<double> A_bd(Ni,Nb);
	A_bd.setFromTriplets(triplets_bd.begin(), triplets_bd.end());
	
	std::cout << "A_int\n" << A_int << std::endl;
	std::cout << "A_bd\n" << A_bd << std::endl;
	
	
	// solving
	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	//solver.isSymmetric(true);	// nice :-)
	solver.compute(A_int);
	Eigen::VectorXd x_int = solver.solve(-A_bd*x_bd);
	Eigen::VectorXd y_int = solver.solve(-A_bd*y_bd);
	
	// back transformation
	xs = x;
	ys = y;
	// interiour nodes
	for(int i = 0; i < Ni; ++i){
		xs(nodes_map[i]) = x_int(i);
		ys(nodes_map[i]) = y_int(i);
	}
	//std::cout << "xs\n" << xs << std::endl;
	//std::cout << "ys\n" << ys << std::endl;
}



// code is ported from Matlab --> not very inteligent...


int main(){
	
	// 
	//	Initalization of mesh
	//
  	Eigen::VectorXd x(10), y(10);
	x << 1.0,0.60,0.12,0.81,0.63,0.09,0.27,0.54,0.95,0.96;
	y << 0.15,0.97,0.95,0.48,0.80,0.14,0.42,0.91,0.79,0.95;
	
	// specify triangles through indices of their vertices
	Eigen::MatrixXi T(11,3);
	T << 7, 1, 2,   5, 6, 2,    4, 1, 7,    6, 7, 2,
		6, 4, 7,   6, 5, 0,    3, 6, 0,    8, 4, 3, 
		3, 4, 6,   8, 1, 4,    9, 1, 8;
/*
	Eigen::MatrixXd T_tmp(T.cast<double>());
	
	mgl::Figure fig1;
	fig1.triplot(T,x,y,"b");
	fig1.save("meshplot_cpp");


	//TriPlot(T_tmp, x, y);
	// Number of nodes of the mesh
	unsigned int N = x.size();
	// Number of triangles of the mesh
	unsigned int M = T.rows();
	
	Eigen::MatrixXi E, Eb, ET;
	processmesh(T,E,Eb);
	getinfo(T, E, ET);
	*/
	Eigen::VectorXd x_ref, y_ref;
	Eigen::MatrixXi T_ref;
	
	refinemesh(x, y, T, x_ref, y_ref, T_ref);
	std::cout << std::endl << "x=\n" << x_ref << std::endl << "y=\n" << y_ref << std::endl <<"T_ref\n" << T_ref << std::endl;
	
	
	//Eigen::MatrixXd T_tmp2(T_ref.cast<double>());
	//TriPlot(T_tmp2, x_ref, y_ref);
	
	Eigen::MatrixXi E, Eb, ET;
	processmesh(T_ref,E,Eb);
	//std::cout << "Eb\n" << Eb << std::endl;
	//std::cout << "E\n" << E << std::endl;
	
	
	mgl::Figure fig3;
	fig3.triplot(T_ref,x_ref,y_ref,"b");
	fig3.save("refine0");
	
	Eigen::VectorXd xs, ys;
	smoothmesh(x_ref, y_ref, T_ref, xs, ys);
	
	mgl::Figure fig2;
	fig2.triplot(T_ref,xs,ys,"b");
	fig2.save("refine1");
	
	
	return 0;
}
