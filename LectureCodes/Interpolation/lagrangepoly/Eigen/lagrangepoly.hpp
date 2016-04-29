# include <Eigen/Dense>

using Eigen::VectorXd;

/* direct evaluation in x of the n^th Lagrange polynomial relative to 'nodes'
 * lagrangepoly(nodes(j), k-1, nodes) = \delta_{j,k} 
 * IN:  x = vector of evaluation points
 *      index n = integer in [0, nodes.size()]
 *      nodes = vector of doubles
 * OUT: y = lagrange polynomial evaluated in x */
void lagrangepoly(const VectorXd& x, const unsigned& n, const VectorXd& nodes, VectorXd& y) {
  Eigen::ArrayXd L = Eigen::ArrayXd::Ones(x.size());
  for (unsigned j = 0; j < nodes.size(); ++j) {
    if (n == j) continue; // avoid division by zero
    L = L * ( x.array() - nodes(j) )/( nodes(n) - nodes(j) );
  }
  y = L.matrix();
}


