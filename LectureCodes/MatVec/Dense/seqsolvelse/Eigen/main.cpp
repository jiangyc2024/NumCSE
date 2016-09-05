#include <iostream>
#include <cmath>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

int main () {
	int N = 10;
	int n = 100;
	VectorXd x(n);
	x.setOnes();
	MatrixXd A(n,n);
	A.setRandom();
	auto some_function = [](VectorXd &x){return x*42;};
	VectorXd b = some_function(x);
	
	
/* SAM_LISTING_BEGIN_0 */
// Setting: \Blue{$N \gg 1$}, 
// large matrix \Blue{$\VA\in\bbK^{n,n}$}
for(int j = 0; j < N; ++j){
	x = A.lu().solve(b);
	b = some_function(x);
}
/* SAM_LISTING_END_0 */


return 0;
}
