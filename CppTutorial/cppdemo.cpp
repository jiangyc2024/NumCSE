// **********************************************************************
// Testing of C++11 festures
// **********************************************************************

#include <cassert>

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;

// Concatenation of sequential containers

template <typename SEQ>
class ConcatIterator {
 public:
  using it_t = typename SEQ::const_iterator;
  using inc_t = typename SEQ::difference_type;
  using si_t = typename std::iterator_traits<it_t>;
  using value_type = typename si_t::value_type;
  using reference = typename si_t::reference;
  using pointer = typename si_t::pointer;
  using difference_type = typename si_t::difference_type;
  using iterator_category = typename si_t::iterator_category;

  ConcatIterator(const SEQ &s1_, const SEQ &s2_, it_t it_)
      : s1(&s1_), s2(&s2_), it(it_) {}
  ConcatIterator(const ConcatIterator &ci) : s1(ci.s1), s2(ci.s2), it(ci.it) {}
  ConcatIterator(ConcatIterator &&) noexcept = default;
  ~ConcatIterator() = default;
  ConcatIterator &operator=(const ConcatIterator &ci) {

    if(this == &ci) { 
      return *this;
    }

    s1 = ci.s1; 
    s2 = ci.s2; 
    it = ci.it;
    return *this;
  }
  ConcatIterator &operator=(ConcatIterator &&) noexcept = default;
  bool operator==(const ConcatIterator<SEQ> &ci) const { return (it == ci.it); }
  bool operator!=(const ConcatIterator<SEQ> &ci) const {
    return !(operator==(ci));
  }
  ConcatIterator &operator++() {
    if (++it == s1->end()) {
      it = s2->begin();
    }
    return *this;
  }
  ConcatIterator &operator++(int) {
    auto tmp = *this;
    operator++();
    return tmp;
  }
  ConcatIterator operator+=(inc_t i) {
    while (i--) {
      operator++();
    }
    return *this;
  }
  value_type operator*() const { return *it; }
  it_t operator->() const { return it; }

 private:

  //std::max_element requires operator=, which cannot be implemented with const &
  const SEQ * s1;
  const SEQ * s2;
  it_t it;
};

template <typename SEQ>
class Concat {
 public:
  using const_iterator = ConcatIterator<SEQ>;
  using value_type = typename SEQ::value_type;
  using reference = typename SEQ::reference;
  using difference_type = typename SEQ::difference_type;
  using size_type = typename SEQ::size_type;
  Concat(const SEQ &s1_, const SEQ &s2_) : s1(s1_), s2(s2_) {}
  [[nodiscard]] const_iterator begin() const {
    return ConcatIterator<SEQ>(s1, s2, s1.begin());
  }
  [[nodiscard]] const_iterator end() const {
    return ConcatIterator<SEQ>(s1, s2, s2.end());
  }
  [[nodiscard]] size_type size() const { return (s1.size() + s2.size()); }

 private:
  const SEQ &s1, &s2;
};

void concattest() {
  using cv_t = Concat<std::vector<double>>;
  std::vector<double> v1 = {1.5, 2.5, 3.5, 4.5};
  std::vector<double> v2 = {7.0, 8.0, 9.0, 10.0};
  cv_t cv(v1, v2);
  for (typename cv_t::const_iterator it = cv.begin(); it != cv.end(); ++it){ //NOLINT(modernize-loop-convert)
    cout << *it << ", " << endl;
  }
  cout << "range loop: ";
  for (auto &&x : cv) {
    cout << x << ", ";
  }
  cout << endl;
  double s = 0.0;
  std::for_each(cv.begin(), cv.end(), [&s](double x) { s += x; });
  cout << "sum = " << s << endl;
  cout << "max = " << *std::max_element(cv.begin(), cv.end()) << endl;
}

// Template deduction through constructor
template <typename T>
class MyClsTempl {
 public:
  using type_t = T;
  MyClsTempl() { ptr = nullptr; }
  explicit MyClsTempl(T &x) { ptr = &x; }
  template <typename U>
  [[nodiscard]] T memfn(const T &x, const U &y) const {
    return x == y ? x : *ptr;
  }

 private:
  T *ptr;
};

// Function template
template <typename ScalarType, typename VectorType>
VectorType saxpy(ScalarType alpha, const VectorType &x, const VectorType &y) {
  return (alpha * x + y);
}

// Function type wrappers

double binop(double arg1, double arg2) { return (arg1 / arg2); }

void functionwrapper() {
  std::vector<std::function<double(double, double)>> fnvec;
  fnvec.emplace_back(binop);
  fnvec.emplace_back([](double x, double y) -> double { return y / x; });

  for (auto const & fn : fnvec) {
    std::cout << fn(3, 2) << std::endl;
  }
}

// Initializer list
void initlist(int n = 42) {
  for (int64_t i : {2, 3, 4, 5, 6, 7, 8, 9, n}) {
    std::cout << i << std::endl;
  }
}

// foreach facility
template <typename T>
void printElements(const T &coll) {
  for (const auto &elem : coll) {
    std::cout << elem << std::endl;
  }
}

//using gsl::owner without actual Guideline Support Library
namespace gsl { 

  template<class T>
  using owner = T; 
} //namespace gsl

// move semantics
class X {
 public:
  X() = default;

  // Range initialization
  template <typename It>
  X(It begin, It end) {
    using value_t = typename It::value_type;
    std::vector<value_t> tmp;
    for (auto it = begin; it != end; ++it, n++) {
      tmp.push_back(*it);
    }
    data = new double[n];
    for (int l = 0; l < n; l++) { 
      data[l] = tmp[l];
    }
  }

  // Initialization from a collection
  template <typename Coll>
  explicit X(const Coll &c) : n(c.size()), data(new double[n]) {
    size_t i = 0;
    for (auto x : c) {
      data[ i ] = x;
      i ++;
    }
  }
  // Copy constructor
  X(const X &x) : n(x.n), data(new double[n]) {
    cout << "X copy constructor: n = " << n << endl;
    for (int l = 0; l < n; l++) {
      data[l] = x.data[l];
    }
  }

  // Move constructor
  X(X &&x) noexcept : n(x.n), data(x.data) {
    cout << "X move constructor: n = " << n << endl;
    x.n = 0;
    x.data = nullptr;
  }

  virtual ~X() {
    cout << "X destructor: n = " << n << endl;
    delete[] data;
  }

  X &operator=(const X &x) {
    cout << "X assignment: n = " << x.n << endl;
    if(this == &x) { 
      return *this;
    }
    delete[] data;
    n = x.n;
    data = new double[n];
    for (int l = 0; l < n; l++) {
      data[l] = x.data[l];
    }
    return *this;
  }

  X &operator=(X &&x) noexcept {
    cout << "X assign & move: n = " << n << ", x.n = " << x.n << endl;
    delete[] data;
    // Stealing the pointers from Rvalue object
    n = x.n;
    data = x.data;
    x.n = 0;
    x.data = nullptr;
    return *this;
  }

  friend std::ostream &operator<<(std::ostream &o, const X &x);

 private:
  int n {0};
  gsl::owner<double*> data {nullptr};
};

std::ostream &operator<<(std::ostream &o, const X &x) {
  o << '[';
  for (int l = 0; l < x.n; l++) {
    o << x.data[l] << ' ';
  }
  return o << "] ";
}

template <typename T>
std::tuple<T, T, std::vector<T>> extcumsum(const std::vector<T> &v) {
  // Local summation variable captured by reference by lambda function
  T sum{};
  // temporary vector for returning cumulative sum
  std::vector<T> w{};
  // cumulative summation
  std::transform(v.cbegin(), v.cend(), back_inserter(w), [&sum](T x) {
    sum += x;
    return (sum);
  });
  return (std::make_tuple(*std::min_element(v.cbegin(), v.cend()),
                          *std::max_element(v.cbegin(), v.cend()),
                          std::move(w)));
}

void extcumsumdemo() {
  // initialize a vector from an initializer list
  std::vector<double> v({1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8});
  // Variables for return values
  double minv(0);
  double maxv(0);       // Extremal elements
  std::vector<double> cs;  // Cumulative sums
  std::tie(minv, maxv, cs) = extcumsum(v);
  cout << "min = " << minv << ", max = " << maxv << endl;
  cout << "cs = [ ";
  for (double x : cs) {
    cout << x << ' ';
  }
  cout << "]" << endl;
}

void newextcumsumdemo() {
  // initialize a vector from an initializer list
  std::vector<double> v({1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8});
  // Definition and assignment of multiple return variables
  auto [minv, maxv, cs] = extcumsum(v);
  cout << "min = " << minv << ", max = " << maxv << endl;
  cout << "cs = [ ";
  for (double x : cs) {
    cout << x << ' ';
  }
  cout << "]" << endl;
}

// Returning several objects
std::tuple<std::string, X, std::size_t> multireturn(const std::vector<double> &v) {
  std::tuple<std::string, X, std::size_t> t("tuple", X(v.begin(), v.end()),
                                       v.size());
  return t;
}

void tupletest_copy() {
  cout << "Testing tuples: copy" << endl;
  std::vector<double> v({1.2, 2.3, 3.4, 4.5});
  // Variables for return values
  std::string s{};
  X x{};
  std::size_t sz(0);

  std::tuple<std::string, X, std::size_t> t(multireturn(v));
  cout << "X(tuple) = " << std::get<1>(t) << endl;
  std::tie(s, x, sz) = t;
  cout << "s = " << s << ", x = " << x << ", sz = " << sz << endl;
}

void tupletest_move() {
  cout << "Testing tuples: move" << endl;
  std::vector<double> v({1.2, 2.3, 3.4, 4.5});
  // Variables for return values
  std::string s{};
  X x{};
  std::size_t sz(0);

  //NOLINTNEXTLINE clang-diagnostic-pessimizing-move
  std::tie(s, x, sz) = std::move(multireturn(v));
  cout << "s = " << s << ", x = " << x << ", sz = " << sz << endl;
}

// Demonstration of lambda function and transform algorithm
void lambdademo() {
  // initialize a vector from an initializer list
  std::vector<double> v({1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8});
  // A vector of the same length
  std::vector<double> w(v.size());
  // Do cumulative summation of v and store result in w
  double sum = 0;
  std::transform(v.begin(), v.end(), w.begin(), [&sum](double x) {
    sum += x;
    return sum;
  });
  cout << "sum = " << sum << ", w = [ ";
  for (auto x : w) {
    cout << x << ' ';
  }
  cout << ']' << endl;
}

int main(int argc, char **argv) {
  cout << "Testing of C++11 features" << endl;
  int64_t code(0);
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " <selection>" << endl; //NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    code = -1L;
  } else {
    const auto sel = std::strtol(argv[1], nullptr, 10); //NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    switch (sel) {
      case 8: {
        extcumsumdemo();
        break;
      }
      case 7: {
        lambdademo();
        break;
      }
      case 1: {
        initlist();
        break;
      }
      case 2: {
        printElements(std::vector<double>{1.2, 2.3, 4.5});
        break;
      }
      case 3: {
        std::vector<double> v({1.2, 2.3, 3.4, 4.5});
        X x(v);
        X y(std::move(x));
        cout << "x = " << x << "y = " << y << endl; //NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved)
        break;
      }
      case 4: {
        tupletest_copy();
        break;
      }
      case 5: {
        tupletest_move();
        break;
      }
      case 6: {
        double x = 3.14;
        float y = 3.14;
        MyClsTempl<int> inst1;
        MyClsTempl<double> inst2(x);
        cout << inst2.memfn(x, y) << endl;
        // MyClsTempl<int> z = saxpy(x,inst1,inst1);
        break;
      }
      case 9: {
        cout << "Concatenation of iterators" << endl;
        concattest();
        break;
      }
      case 10: {
        functionwrapper();
        break;
      }
      default: {
        cerr << "Invalid selection" << endl;
        exit(-1L);
      }
    }
  }
  return code;
}
