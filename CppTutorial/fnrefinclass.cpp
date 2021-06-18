/***********************************************************************
 *                                                                     *
 * Demo code                                                           *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: R.H.                                                        *
 * Date: April 2021                                                    *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

// Header for basic IO
#include <functional>
#include <iostream>
#include <vector>

namespace FnRefInClass {
template <typename Functor>
class ClassWithFnRef {
 public:
  explicit ClassWithFnRef(Functor& fn) : fn_(fn) {}
  ClassWithFnRef() = delete;
  ClassWithFnRef(const ClassWithFnRef&) = default;
  ClassWithFnRef(ClassWithFnRef&&) = delete;
  ClassWithFnRef& operator=(const ClassWithFnRef&) = default;
  ClassWithFnRef& operator=(ClassWithFnRef&&) = delete;
  ~ClassWithFnRef() = default;

  double operator()(double x) const {
    const double res = fn_(x);
    std::cout << typeid(*this).name() << "(" << x << ") = " << res << std::endl;
    return res;
  }

 private:
  Functor& fn_;
};

auto createFnRefInClass(double a) {
  // A local lambda function
  auto fn = [a](double x) { return a + x; };
  // Create an object
  return ClassWithFnRef<decltype(fn)>(fn);
}

template <typename Functor>
double testPassRefInClass(const ClassWithFnRef<Functor>& obj) {
  // Evaluations
  double s = 0.0;
  for (double x : std::vector<double>({1.0, 1.5, 2.0, 2.5, 3.0})) {
    s += obj(x);
  }
  return s;
}

double testFnRefInClass() {
  // A local variable
  double a = 3.14;
  // Create an object
  auto obj{createFnRefInClass(a)};
  return testPassRefInClass(obj);
}

double testLocalFnRefInClass() {
  // A local variable
  double e = 2.718;
  // A local lambda function
  auto fn = [e](double x) { return e + x; };
  // Create object
  ClassWithFnRef<decltype(fn)> obj(fn);
  // Evaluations
  double s = 0.0;
  for (double x : std::vector<double>({1.0, 1.5, 2.0, 2.5, 3.0})) {
    s += obj(x);
  }
  return s;
}

}  // namespace FnRefInClass

int main(int /*argc*/, char** /*argv*/) {
  std::cout << "OK: A reference to a lambda function is used while still alive"
            << std::endl;
  (void)FnRefInClass::testLocalFnRefInClass();
  std::cout << "\nDANGER: A reference to a lambda function is used out of scope!"
            << std::endl;
  double sum = FnRefInClass::testFnRefInClass();
  std::cout << "sum = " << sum << std::endl;
  return 0;
}
