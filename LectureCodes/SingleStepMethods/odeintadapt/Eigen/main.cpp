#include "figure.hpp"
#include "odeintadapt.hpp"
#include <iostream>
#include <math.h>

int main()
{
  using State_t = double;
  using DiscEvolOp = std::function<State_t(double,State_t)>;
  
  // Differential equation $\dot{y} = y^2$
  auto f = [](double x){ return pow(x,2); };
  auto norm = [](double x){ return fabs(x); };
  
  // explicit euler (order 1) 
  DiscEvolOp psilow = [&](double h, State_t y){ return  y + h*f(y); };
  
  // explicit trapezoidal (order 2)
  DiscEvolOp psihigh = [&](double h, State_t y) {
    double k1 = f(y);
    double k2 = f(y + h*k1);
    return y + (h/2.)*(k1+k2);
  };
  
  
  double y0 = 0.5; 
  std::vector<std::pair<double, double>>
    states = odeintadapt(psilow,psihigh,y0, 1.9, 0.2, 1e-2, 1e-2, 1e-4,norm);

  // Save result for graphical output
  Eigen::VectorXd y(states.size());
  Eigen::VectorXd t(states.size());
  for (int i=0; i<states.size(); ++i)
    {
      t(i) = states[i].first;
      y(i) = states[i].second;
    }

  // Graphical output
  mgl::Figure lin;
  lin.plot(t, y, "#r^").label("y_k");
  lin.fplot("1/(2-x)").label("y(t)").style("b-");
  
  std::cout << "Last state at t="<< states.back().first << " with y=" << states.back().second << std::endl;
  
  // lin.addlabel("y' = y^2", "r");
  lin.legend(1,1);
  lin.xlabel("t");
  lin.ylabel("y");
  lin.title("Simple local stepsize control");
  lin.setFontSize(3);
  lin.save("odeintadapt");
}
