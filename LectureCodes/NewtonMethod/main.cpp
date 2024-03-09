// **********************************************************************
// Eigen codes related to Newton's method
// **********************************************************************

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <span>

template <typename FuncType,typename DervType,typename Scalar>
Scalar newton1D(const FuncType &F,const DervType &DF,
		const Scalar &x0,double rtol,double atol)
{
  Scalar s;
  Scalar z = x0;
  do {
    s = F(z)/DF(z);  // compute Newton correction
    z -= s;           // compute next iterate
  }
  // correction based termination (relative and absolute)
  while ((std::abs(s) > rtol*std::abs(z)) && (std::abs(s) > atol));
  return (z);
}

template <typename FuncType,typename DervType,typename VecType>
void newton(const FuncType &F,const DervType &DFinv,
	    VecType &x,double rtol,double atol)
{
  using index_t = typename VecType::Index;
  using scalar_t = typename VecType::Scalar;
  const index_t n = x.size();
  VecType s(n);
  scalar_t sn;
  do {
    s = DFinv(x,F(x));  // compute Newton correction
    x -= s;                  // compute next iterate
    sn = s.norm();
  }
  // simple correction based termination (relative and absolute)
  while ((sn > rtol*x.norm()) && (sn > atol));
}

template <typename FuncType,typename JacType,typename VecType>
void newton_stc(const FuncType &F,const JacType &DF,
	    VecType &x,double rtol,double atol)
{
  using scalar_t = typename VecType::Scalar;
  scalar_t sn;

  do {
    auto jacfac = DF(x).lu(); // LU-factorize Jacobian
    x -= jacfac.solve(F(x));  // Compute next iterate
    // Compute norm of simplified Newton correction
    sn = jacfac.solve(F(x)).norm();  
    std::cout << "|ss| = " << sn << std::endl;
  }
  // Termination based on simplified Newton correction
  while ((sn > rtol*x.norm()) && (sn > atol));
}

template <typename FuncType,typename JacType,typename VecType>
void dampnewton(const FuncType &F,const JacType &DF,
	        VecType &x,double rtol,double atol)
{
  using index_t = typename VecType::Index;
  using scalar_t = typename VecType::Scalar;
  const index_t n = x.size();
  const scalar_t lmin = 1E-3; // Minimal damping factor
  scalar_t lambda = 1.0;      // Actual damping factor
  VecType s(n);
  VecType st(n);      // Newton corrections
  VecType xn(n);      // Tentative new iterate
  scalar_t sn;
  scalar_t stn;       // Norms of Newton corrections

  do {
    auto jacfac = DF(x).lu(); // LU-factorize Jacobian
    s = jacfac.solve(F(x));   // Newton correction
    sn = s.norm();            // Norm of Newton correction
    lambda *= 2.0; 
    do {
      lambda /= 2;
      if (lambda < lmin) {
        throw std::runtime_error("No convergence: lambda -> 0");
      }
      xn = x-lambda*s;           // Tentative next iterate
      st = jacfac.solve(F(xn));  // Simplified Newton correction  
      stn = st.norm(); 
      std::cout << "Inner: |stn| = " << stn << std::endl;
    }
    while (stn > (1-lambda/2)*sn); // Natural monotonicity test
    x = xn; // Now: xn accepted as new iterate
    lambda = std::min(2.0*lambda,1.0); // Try to mitigate damping
    std::cout << "|sn| = " << sn << ", x = [" << x.transpose() << "]" << std::endl;
  }
  // Termination based on simplified Newton correction
  while ((stn > rtol*x.norm()) && (stn > atol));
}




void newton1Ddriver(double x0)
{
  auto F = [](double x) { return x*exp(x)-1.0; };
  auto DF = [](double x) { return exp(x)*(1.0+x); };
  const double z = newton1D(F,DF,x0,1E-6,1E-8);
  std::cout << "F(z) = " << F(z) << std::endl;
}

void newton2Ddriver()
{
  auto F = [](const Eigen::Vector2d &x) {
    Eigen::Vector2d z;
    const double x1 = x(0);
    const double x2 = x(1);
    z << x1*x1-2*x1-x2+1,x1*x1+x2*x2-1;
    return(z);
  };
  auto DF = [](const Eigen::Vector2d &x,const Eigen::Vector2d &f) {
    Eigen::Matrix2d J;
    const double x1 = x(0);
    const double x2 = x(1);
    J << 2*x1-2,-1,2*x1,2*x2;
    Eigen::Vector2d s = J.lu().solve(f);
    return(s);
  };
  Eigen::Vector2d x; x << 2,3;
  newton(F,DF,x,1E-6,1E-8);
  std::cout << "||F(x)|| = " << F(x).norm() << std::endl;
}

void newtonstc2Ddriver()
{
  auto F = [](const Eigen::Vector2d &x) {
    Eigen::Vector2d z;
    const double x1 = x(0);
    const double x2 = x(1);
    z << x1*x1-2*x1-x2+1,x1*x1+x2*x2-1;
    return(z);
  };
  auto Jac = [](const Eigen::Vector2d &x) {
    Eigen::Matrix2d J;
    const double x1 = x(0);
    const double x2 = x(1);
    J << 2*x1-2,-1,2*x1,2*x2;
    return(J);
  };
  Eigen::Vector2d x; x << 2,3;
  newton_stc(F,Jac,x,1E-6,1E-8);
  std::cout << "||F(x)|| = " << F(x).norm() << std::endl;
}

void dampnewtondriver(double x0 = 2,double x1 =3)
{
  auto F = [](const Eigen::Vector2d &x) {
    Eigen::Vector2d z;
    const double x1 = x(0);
    const double x2 = x(1);
    z << x1*x1-2*x1-x2+1,x1*x1+x2*x2-1;
    return(z);
  };
  auto Jac = [](const Eigen::Vector2d &x) {
    Eigen::Matrix2d J;
    const double x1 = x(0);
    const double x2 = x(1);
    J << 2*x1-2,-1,2*x1,2*x2;
    return(J);
  };
  Eigen::Vector2d x; x << x0,x1;
  std::cout << "DAMPED NEWTON, initial guess x = ["
	    << x.transpose() << "]" << std::endl;
  try { dampnewton(F,Jac,x,1E-6,1E-8); }
  catch(char const *msg) { std::cout << msg << std::endl; }
  std::cout << "||F(x)|| = " << F(x).norm() << std::endl;
}

//NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc,char **argv)
{
  auto args = std::span(argv, argc);
  int exitCode = 0;
  std::cout << "EIGEN NEWTON METHOD codes" << std::endl;
  
  if(argc > 1){
    const int64_t sel = std::strtol(args[1], nullptr, 10);
    switch (sel) {
    case 1: { newton1Ddriver(2.0); break; }
    case 2: { newton2Ddriver(); break; } 
    case 3: { newtonstc2Ddriver(); break; } 
    case 4: {
      double x0 = 2;
      double x1 = 3;
      if (argc == 4) {
        x0 = std::strtod(args[2], nullptr);
        x1 = std::strtod(args[3], nullptr);
      }
      dampnewtondriver(x0,x1);
      break; 
    } 
    default: { std::cerr << "Invalid selection" << std::endl; exitCode = -1L; }
    }
  }
  else {
    std::cerr << "Usage: " << args[0] << " <selection>" << std::endl;
    exitCode = -1L;
  }
  return exitCode;
}
