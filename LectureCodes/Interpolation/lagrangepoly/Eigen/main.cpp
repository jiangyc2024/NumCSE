# include "./lagrangepoly.hpp"
# include <figure/figure.hpp>

int main () {
  Eigen::VectorXd nodes(4);
  nodes << -0.8, 0, 0.3, 1;

  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(200, -1, 1),
                  l0, l1, l2, l3;

  lagrangepoly(x, 0, nodes, l0);
  lagrangepoly(x, 1, nodes, l1);
  lagrangepoly(x, 2, nodes, l2);
  lagrangepoly(x, 3, nodes, l3);

  mgl::Figure fig;
  fig.plot(nodes, Eigen::VectorXd::Zero(nodes.size()), " ko").label("Nodes");
  fig.plot(x, l0).label("L_0(x)");
  fig.plot(x, l1).label("L_1(x)");
  fig.plot(x, l2).label("L_2(x)");
  fig.plot(x, l3).label("L_3(x)");
  fig.legend();
  fig.save("poly");

  return 0;
}
