///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <vector>
#include <Eigen/Dense>
#include <iostream>

// explicit ode solver based on the Dormand–Prince method and predictive stepsize control
template <class Func, class State, class NormFunc>
std::vector<std::pair<double, State>> rkf45(Func &f, State &initialState, double T,
					    NormFunc &norm, double h0 = 1e-3, double reltol= 1e-7, double abstol= 1e-7,
					    double hmin= 1e-4, bool predictiveStepsize = true)
{
  double h = h0; // initial stepsize
  double t = 0;  // initial time 
  
  std::vector<std::pair<double, State>> states; // saves the state vector after every steps
  states.push_back({t, initialState});

  // Definition of the embedded Runge-Kutta method: Butcher tableau
  Eigen::VectorXd E(7);
  Eigen::MatrixXd A(7,7);
  E << 5179./57600.,0,7571./16695.,393./640.,-92097./339200.,187./2100.,1./40.;
  A << 0,0,0,0,0,0,0,
    1./5.,0,0,0,0,0,0,
    3./40.,9./40., 0,0,0,0,0,
    44./45.,-56./15., 32./9.,0,0,0,0,
    19372./6561.,-25360./2187.,64448./6561., -212./729.,0,0,0,
    9017./3168.,-355/33,46732./5247., 49./176.,-5103./18656.,0,0,
    35./384., 0, 500./1113.,125./192.,-2187./6784.,11./84.,0;
	
  State y0, y1, y1low;
  y0 = initialState;

  while ((states.back().first < T) && (h >= hmin)) {
    y1 = y0;
    y1low = y0;
    std::vector<State> K; K.reserve(7);

    for (unsigned int i=0; i<7; ++i) {
      // \Blue{$\Vy_{0}+\tstep \sum\limits_{j=1}^{i-1}a_{ij}\Vk_{j})$}
      State incr = y0;
      for (unsigned int j=0; j<i; ++j) {
	incr += h*A(i,j)*K.at(j);
      }
      K.push_back(f(incr));
      // \Blue{$\Vy_{1}  += \tstep b_{i}\Vk_{i}$}
      y1 += h*A(6,i)*K.back(); // fifth-order solution
      y1low += h*E(i)*K.back(); // fourth-order solution
    }
    // stepsize control
    double est = norm(y1-y1low); // estimated error
    if (predictiveStepsize)
      {
	double tol = std::max(reltol*norm(y0), abstol);
	double hOld = h;
	// Optimal stepsize according to \eqref{eq:ssc}\Label[line]{ssctrl:6b} with \Blue{$p=4$}
	h = h*std::max(0.5, std::min(2., pow(tol/est, 1./5.))); 
	if (est < tol || h == hOld) {
	  //step \Magenta{accepted}
	  y0 = y1;
	  states.push_back({states.back().first+std::min(T-states.back().first, h), y0});
	  if (h==hOld && est >= tol)
	    std::cerr << "Warning: Unable to meet tolerance requirements at t=" << states.back().first << std::endl;
	}
      }	
    else
      {
	if (est < std::max(reltol*norm(y0), abstol)) // \label{odeintadapt:6}
	  {
	    y0 = y1; 
	    states.push_back({states.back().first +	std::min(T-states.back().first, h), y0}); 
	    h = 1.1*h;  // step \Magenta{accepted}, try with increased stepsize 
	  }
	else
	  {
	    h = h/2; // step \Magenta{rejected}, try with half the stepsize 
	  }
      }
  }
  
  if (h < hmin)
    {
      std::cerr 	<< "Warning: Failure at t=" 
			<< states.back().first 
			<< ". Unable to meet integration tolerances without reducing the step size below the smallest value allowed ("<< hmin <<") at time t." << std::endl;
    }
  
  return states;
  
}
