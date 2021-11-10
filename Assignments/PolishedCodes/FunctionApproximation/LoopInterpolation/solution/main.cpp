#include "adaptiveintp.hpp"

#define PI M_PI

int main() {
  auto c = [](double t) -> Vector2d {
    Vector2d ct;
    ct << std::cos(2 * PI * t) + 2. / 3. * std::cos(4 * PI * t),
        3. / 2. * std::sin(2 * PI * t);
    return ct;
  };

  VectorXd t = VectorXd::LinSpaced(4, 0, 0.75);

  std::vector<Vector2d> Sigma(t.size());
  for (unsigned int i = 0; i < t.size(); ++i) {
    Sigma[i] = c(t(i));
  }

  unsigned int neval = 5;
  VectorXd x(neval);
  x << 0.0, 0.4, 0.7, 0.9, 1.;
  /*
   *  run closedPolygonalInterpolant
   */
  std::vector<Vector2d> v = closedPolygonalInterpolant(Sigma, x);
  std::cout << "\nResult of closedPolygonalInterpolant:\n";
  for (auto & i : v) {
    std::cout << i.transpose() << std::endl;
  }

  /*
   *  test of closedPolygonalInterpolant
   */
  std::vector<Vector2d> vtest;

  vtest.emplace_back(1.66667, 0);
  vtest.emplace_back(-0.520348, 0.841566);
  vtest.emplace_back(-0.508862, -1.39855);
  vtest.emplace_back(0.941491, -0.466185);
  vtest.emplace_back(1.66667, 0);

  if (vtest.size() == v.size()) {
    double err = 0;
    for (unsigned int i = 0; i < v.size(); ++i) {
      err += (vtest[i] - v[i]).norm();
    }
    if (err < 1e-4)
      std::cout << "Test for closedPolygonalInterpolant passed!\n\n";
    else
      std::cout
          << "Test for closedPolygonalInterpolant failed: wrong output.\n\n";
  } else
    std::cout << "Test for closedPolygonalInterpolant failed: wrong size.\n\n";

  /*
   *  run closedHermiteInterpolant
   */
  std::vector<Vector2d> w = closedHermiteInterpolant(Sigma, x);
  std::cout << "\nResult of closedHermiteInterpolant:\n";
  for (auto & i : w) {
    std::cout << i.transpose() << std::endl;
  }

  /*
   *  test of closedHermiteInterpolant
   */
  std::vector<Vector2d> wtest;

  wtest.emplace_back(1.66667, 0);
  wtest.emplace_back(-0.603704, 0.853533);
  wtest.emplace_back(-0.579763, -1.64522);
  wtest.emplace_back(1.19441, -0.927398);
  wtest.emplace_back(1.66667, 0);

  if (wtest.size() == w.size()) {
    double err = 0;
    for (unsigned int i = 0; i < w.size(); ++i) {
      err += (wtest[i] - w[i]).norm();
    }
    if (err < 1e-4)
      std::cout << "Test for closedHermiteInterpolant passed!\n\n";
    else
      std::cout
          << "Test for closedHermiteInterpolant failed: wrong output.\n\n";
  } else
    std::cout << "Test for closedHermiteInterpolant failed: wrong size.\n\n";

  /*
   *  run adaptedHermiteInterpolant
   */

  std::pair<std::vector<Vector2d>, std::vector<Vector2d> > p =
      adaptedHermiteInterpolant(c, 5, x, 1e-3);
  std::cout << "\nResult of adaptedHermiteInterpolant:\n";
  for (auto & i : p.first) {
    std::cout << i.transpose() << std::endl;
  }

  /*
   *  test of adaptedHermiteInterpolant
   */

  std::vector<Vector2d> wwtest;

  wwtest.emplace_back(1.66667, 0);
  wwtest.emplace_back(-0.606359, 0.88769);
  wwtest.emplace_back(-0.536755, -1.49042);
  wwtest.emplace_back(1.17946, -0.765441);
  wwtest.emplace_back(1.66667, -3.67761e-16);

  if (wwtest.size() == (p.first).size()) {
    double err = 0;
    for (unsigned int i = 0; i < wwtest.size(); ++i) {
      err = std::max(err, (wwtest[i] - p.first[i]).norm());
    }
    if (err < 1e-3)
      std::cout << "Test for adaptedHermiteInterpolant passed!\n\n";
    else
      std::cout
          << "Test for adaptedHermiteInterpolant failed: wrong output.\n\n";
  } else
    std::cout << "Test for adaptedHermiteInterpolant failed: wrong size.\n\n";
  /*
   *  test of plotKite
   */
  plotKite("./cx_out/fig1.png", 1.0e-1);
  plotKite("./cx_out/fig2.png", 1.0e-2);
  plotKite("./cx_out/fig3.png", 1.0e-3);
}
