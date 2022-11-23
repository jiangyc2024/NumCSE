// **********************************************************************
// Simple custom vector class
// **********************************************************************

#include <algorithm>
#include <list>

// **********************************************************************
// Simple custom vector class
// **********************************************************************

#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

using std::cout;
using std::endl;

namespace myvec {
class MyVector {
public:
  using value_t = double;
  // Constructor creating constant vector, also default constructor
  explicit MyVector(std::size_t n = 0, double val = 0.0);
  // Constructor: initialization from an STL container
  template <typename Container> MyVector(const Container &v);
  // Constructor: initialization from an STL iterator range
  template <typename Iterator> MyVector(Iterator first, Iterator last);
  // Copy constructor, computational cost O(n)
  MyVector(const MyVector &mv);
  // Move constructor, computational cost O(1)
  MyVector(MyVector &&mv);
  // Assignment operator, computational cost O(n)
  MyVector &operator=(const MyVector &mv);
  // Move assignment operator, computational cost O(1)
  MyVector &operator=(MyVector &&mv);
  // Destructor
  virtual ~MyVector(void);
  // Type conversion to STL vector
  operator std::vector<double>() const;

  // Returns length of vector
  std::size_t size(void) const { return n; }
  // Access operators: rvalue \& lvalue, with range check
  double operator[](std::size_t i) const;
  double &operator[](std::size_t i);
  // Comparison operators
  bool operator==(const MyVector &mv) const;
  bool operator!=(const MyVector &mv) const;
  // Transformation of a vector by a function of signature double foo(double)
  template <typename Functor> MyVector &transform(Functor &&f);

  // Overloaded arithmetic operations
  // In place vector addition: x += y;
  MyVector &operator+=(const MyVector &mv);
  // In place vector subtraction: x-= y;
  MyVector &operator-=(const MyVector &mv);
  // In place scalar multiplication: x *= a;
  MyVector &operator*=(double alpha);
  // In place scalar division: x /= a;
  MyVector &operator/=(double alpha);
  // Vector addition
  MyVector operator+(MyVector mv) const;
  // Vector subtraction
  MyVector operator-(const MyVector &mv) const;
  // Scalar multiplication from right and left: x = a*y; x = y*a
  MyVector operator*(double alpha) const;
  friend MyVector operator*(double alpha, const MyVector &);
  // Scalar divsion: x = y/a;
  MyVector operator/(double alpha) const;
  // Euclidean norm
  double norm(void) const;
  // Euclidean inner product
  double operator*(const MyVector &)const;
  // Output function
  friend std::ostream &operator<<(std::ostream &, const MyVector &mv);

  static bool dbg; // Flag for verbose output
private:
  std::size_t n; // Length of vector
  double *data;  // data array (standard C array)
};
} // namespace myvec

// **********************************************************************
// Simple custom vector class
// **********************************************************************

namespace myvec {
bool MyVector::dbg = true;

// Implementation of member functions
MyVector::MyVector(std::size_t _n, double _a) : n(_n), data(nullptr) {
  if (dbg)
    cout << "{Constructor MyVector(" << _n << ") called" << '}' << endl;
  if (n > 0)
    data = new double[_n];
  for (std::size_t l = 0; l < n; ++l)
    data[l] = _a;
}

template <typename Container>
MyVector::MyVector(const Container &v) : n(v.size()), data(nullptr) {
  if (dbg)
    cout << "{MyVector(length " << n << ") constructed from container" << '}'
         << endl;
  if (n > 0) {
    double *tmp = (data = new double[n]);
    for (auto i : v)
      *tmp++ = i;
  }
}

template <typename Iterator>
MyVector::MyVector(Iterator first, Iterator last) : n(0), data(nullptr) {
  n = std::distance(first, last);
  if (dbg)
    cout << "{MyVector(length " << n << ") constructed from range" << '}'
         << endl;
  if (n > 0) {
    data = new double[n];
    std::copy(first, last, data);
  }
}

MyVector::MyVector(const MyVector &mv) : n(mv.n), data(nullptr) {
  if (dbg)
    cout << "{Copy construction of MyVector(length " << n << ")" << '}' << endl;
  if (n > 0) {
    data = new double[n];
    std::copy_n(mv.data, n, data);
  }
}

MyVector::MyVector(MyVector &&mv) : n(mv.n), data(mv.data) {
  if (dbg)
    cout << "{Move construction of MyVector(length " << n << ")" << '}' << endl;
  mv.data = nullptr;
  mv.n = 0;
}

MyVector &MyVector::operator=(const MyVector &mv) {
  if (dbg)
    cout << "{Copy assignment of MyVector(length " << n << "<-" << mv.n << ")"
         << '}' << endl;
  if (this == &mv)
    return (*this);
  if (n != mv.n) {
    n = mv.n;
    if (data != nullptr)
      delete[] data;
    if (n > 0)
      data = new double[n];
    else
      data = nullptr;
  }
  if (n > 0)
    std::copy_n(mv.data, n, data);
  return (*this);
}

MyVector &MyVector::operator=(MyVector &&mv) {
  if (dbg)
    cout << "{Move assignment of MyVector(length " << n << "<-" << mv.n << ")"
         << '}' << endl;
  if (data != nullptr)
    delete[] data;
  n = mv.n;
  data = mv.data;
  mv.n = 0;
  mv.data = nullptr;
  return (*this);
}

MyVector::~MyVector(void) {
  if (dbg)
    cout << "{Destructor for MyVector(length = " << n << ")" << '}' << endl;
  if (data != nullptr)
    delete[] data;
}

MyVector::operator std::vector<double>() const {
  if (dbg)
    cout << "{Conversion to std::vector, length = " << n << '}' << endl;
  return (std::vector<double>(data, data + n));
}

double MyVector::operator[](std::size_t i) const {
  if (i >= n)
    throw(std::logic_error("[] out of range"));
  return data[i];
}

double &MyVector::operator[](std::size_t i) {
  if (i >= n)
    throw(std::logic_error("[] out of range"));
  return data[i];
}

bool MyVector::operator==(const MyVector &mv) const {
  if (dbg)
    cout << "{Comparison ==: " << n << " <-> " << mv.n << '}' << endl;
  if (n != mv.n)
    return (false);
  else {
    for (std::size_t l = 0; l < n; ++l)
      if (data[l] != mv.data[l])
        return (false);
  }
  return (true);
}

bool MyVector::operator!=(const MyVector &mv) const { return !(*this == mv); }

template <typename Functor> MyVector &MyVector::transform(Functor &&f) {
  for (std::size_t l = 0; l < n; ++l)
    data[l] = f(data[l]);
  return (*this);
}

MyVector &MyVector::operator+=(const MyVector &mv) {
  if (dbg)
    cout << "{operator +=, MyVector of length " << n << '}' << endl;
  if (n != mv.n)
    throw(std::logic_error("+=: vector size mismatch"));
  for (std::size_t l = 0; l < n; ++l)
    data[l] += mv.data[l];
  return (*this);
}

MyVector &MyVector::operator-=(const MyVector &mv) {
  if (dbg)
    cout << "{operator -=, MyVector of length " << n << '}' << endl;
  if (n != mv.n)
    throw(std::logic_error("-=: vector size mismatch"));
  for (std::size_t l = 0; l < n; ++l)
    data[l] -= mv.data[l];
  return (*this);
}

MyVector &MyVector::operator*=(double alpha) {
  if (dbg)
    cout << "{operator *=, MyVector of length " << n << '}' << endl;
  for (std::size_t l = 0; l < n; ++l)
    data[l] *= alpha;
  return (*this);
}

MyVector &MyVector::operator/=(double alpha) {
  if (dbg)
    cout << "{operator *=, MyVector of length " << n << '}' << endl;
  for (std::size_t l = 0; l < n; ++l)
    data[l] /= alpha;
  return (*this);
}

MyVector MyVector::operator+(MyVector mv) const {
  if (dbg)
    cout << "{operator +, MyVector of length " << n << '}' << endl;
  if (n != mv.n)
    throw(std::logic_error("+: vector size mismatch"));
  mv += *this;
  return (mv);
}

MyVector MyVector::operator-(const MyVector &mv) const {
  if (dbg)
    cout << "{operator +, MyVector of length " << n << '}' << endl;
  if (n != mv.n)
    throw(std::logic_error("+: vector size mismatch"));
  MyVector tmp(*this);
  tmp -= mv;
  return (tmp);
}

MyVector MyVector::operator*(double alpha) const {
  if (dbg)
    cout << "{operator *a, MyVector of length " << n << '}' << endl;
  MyVector tmp(*this);
  tmp *= alpha;
  return (tmp);
}

MyVector operator*(double alpha, const MyVector &mv) {
  if (MyVector::dbg)
    cout << "{operator a*, MyVector of length " << mv.n << '}' << endl;
  MyVector tmp(mv);
  tmp *= alpha;
  return (tmp);
}

MyVector MyVector::operator/(double alpha) const {
  if (dbg)
    cout << "{operator /, MyVector of length " << n << '}' << endl;
  MyVector tmp(*this);
  tmp /= alpha;
  return (tmp);
}

double MyVector::norm(void) const {
  if (dbg)
    cout << "{norm: MyVector of length " << n << '}' << endl;
  double s = 0;
  for (std::size_t l = 0; l < n; ++l)
    s += (data[l] * data[l]);
  return (std::sqrt(s));
}

double MyVector::operator*(const MyVector &mv) const {
  if (dbg)
    cout << "{dot *, MyVector of length " << n << '}' << endl;
  if (n != mv.n)
    throw(std::logic_error("dot: vector size mismatch"));
  double s = 0;
  for (std::size_t l = 0; l < n; ++l)
    s += (data[l] * mv.data[l]);
  return (s);
}

std::ostream &operator<<(std::ostream &o, const MyVector &mv) {
  o << "[ ";
  for (std::size_t l = 0; l < mv.n; ++l)
    o << mv.data[l] << (l == mv.n - 1 ? ' ' : ',');
  return (o << "]");
}
} // namespace myvec

using namespace myvec;

template <typename Vec>
std::vector<Vec> gramschmidt(const std::vector<Vec> &A, double eps = 1E-14) {
  const int k = A.size();    // no. of vectors to orthogonalize
  const int n = A[0].size(); // length of vectors
  cout << "gramschmidt orthogonalization for " << k << ' ' << n << "-vectors"
       << endl;
  std::vector<Vec> Q({A[0] / A[0].norm()}); // output vectors
  for (int j = 1; (j < k) && (j < n); ++j) {
    Q.push_back(A[j]);
    for (int l = 0; l < j; ++l)
      Q.back() -= (A[j] * Q[l]) * Q[l];
    if (Q.back().norm() < eps * A[j].norm()) { // premature termination ?
      Q.pop_back();
      break;
    }
    Q.back() /= Q.back().norm(); // normalization
  }
  return (Q); // return at end of local scope
}

// Initialization of a sequence of vectors
template <typename Functor>
std::vector<MyVector> initvectors(std::size_t n, std::size_t k, Functor &&f) {
  std::vector<MyVector> A{};
  for (int j = 0; j < (int)k; ++j) {
    A.push_back(MyVector(n));
    for (int i = 0; i < (int)n; ++i) {
      (A.back())[i] = f(i, j);
    }
  }
  return (A);
}

struct SimpleFunction {
  SimpleFunction(double _a = 1.0) : cnt(0), a(_a) {}
  double operator()(double x) {
    cnt++;
    return (x + a);
  }
  int cnt;        // internal counter
  const double a; // increment value
};

int main(int argc, char **argv) {
  cout << "MyVector class implementation" << endl;
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <selection>" << endl;
    std::cerr << "1 : plain allocatiobn" << std::endl;
    std::cerr << "2 : initialization from list" << std::endl;
    std::cerr << "3 : initialization from container" << std::endl;
    std::cerr << "4 : entrywise operation" << std::endl;
    std::cerr << "5 : Gram-Schmidt" << std::endl;
    return (-1L);
  } else {
    try {
      const int sel = atoi(argv[1]);
      switch (sel) {
      case 1: {
        MyVector mv(std::size_t(10));
        cout << "Zero vector: " << mv << endl;
        break;
      }
      case 2: {
        std::list<int> lst = {1, 2, 3, 5, 7, 11, 13};
        MyVector mv(lst);
        cout << "Initialized from list: " << mv << endl;
        cout << "mv[5] = " << mv[5] << endl;
        mv[1] = 3.14;
        cout << "Modified: " << mv << endl;
        cout << "Illegal: " << mv[17] << endl;
        break;
      }
      case 3: {
        myvec::MyVector::dbg = true;
        std::vector<int> ivec = {1, 2, 3, 5, 7, 11, 13};
        myvec::MyVector v1(ivec.cbegin(), ivec.cend());
        myvec::MyVector v2(ivec);
        myvec::MyVector vr(ivec.crbegin(), ivec.crend());
        cout << "v1 = " << v1 << endl;
        cout << "v2 = " << v2 << endl;
        cout << "vr = " << vr << endl;
        break;
      }
      case 6: {
        myvec::MyVector::dbg = true;
        myvec::MyVector v1(
            std::vector<double>({1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9}));
        myvec::MyVector v2(2.0 * v1);
        myvec::MyVector v3(std::move(v1));
        cout << "v1 = " << v1 << endl;
        cout << "v2 = " << v2 << endl;
        cout << "v3 = " << v3 << endl;
        break;
      }
      case 7: {
        myvec::MyVector::dbg = true;
        myvec::MyVector x(
            std::vector<double>({1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9}));
        myvec::MyVector y(
            std::vector<double>({2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8}));
        auto z = x + (x * y) * x + 2.0 * y / (x - y).norm();
        break;
      }
      case 4: {
        myvec::MyVector::dbg = false;
        double a = 2.0; // increment
        int cnt = 0;    // external counter used by lambda function
        myvec::MyVector mv(
            std::vector<double>({1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9}));
        mv.transform([a, &cnt](double x) {
          cnt++;
          return (x + a);
        });
        cout << cnt << " operations, mv transformed = " << mv << endl;
        SimpleFunction trf(a);
        mv.transform(trf);
        cout << trf.cnt << " operations, mv transformed = " << mv << endl;
        mv.transform(SimpleFunction(-4.0));
        cout << "Final vector = " << mv << endl;
        break;
      }
      case 5: {
        MyVector::dbg = false;
        const int n = 7;
        const int k = 7;
        std::vector<myvec::MyVector> A(initvectors(
            n, k, [](int i, int j) { return std::min(i + 1, j + 1); }));
        std::vector<myvec::MyVector> Q(gramschmidt(A));
        cout << "Set of vectors to be orthonormalized:" << endl;
        for (const auto &a : A) {
          cout << a << endl;
        }
        cout << "Output of Gram-Schmidt orthonormalization: " << endl;
        for (const auto &q : Q) {
          cout << q << endl;
        }
        cout << "Testing orthogonality:" << endl;
        for (const auto &qi : Q) {
          for (const auto &qj : Q)
            cout << std::setprecision(3) << std::setw(9) << qi * qj << ' ';
          cout << endl;
        }
        break;
      }

      default: {
        std::cerr << "Invalid selection" << endl;
        exit(-1L);
      }
      }
    } catch (std::exception &e) {
      std::cerr << "ERROR: " << e.what() << endl;
    }
  }
  return (0);
}
