// **********************************************************************
// Testing of C++11 festures
// **********************************************************************

#include <assert.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <utility>
#include <functional>

using namespace std;

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
  
  ConcatIterator(const SEQ &s1_,const SEQ &s2_,it_t it_):s1(s1_),s2(s2_),it(it_) {}
  ConcatIterator(const ConcatIterator &ci):s1(ci.s1),s2(ci.s2),it(ci.it) {}
  ConcatIterator &operator = (const ConcatIterator &ci) {
    // This does not compile
    // s1=ci.s1; s2=ci.s2; it = ci.it;
    return *this; }
  bool operator == (const ConcatIterator<SEQ> &ci) const { return (it == ci.it); }
  bool operator != (const ConcatIterator<SEQ> &ci) const { return !(operator==(ci)); }
  ConcatIterator &operator ++ () { if (++it == s1.end()) it = s2.begin(); return *this; }
  ConcatIterator &operator ++ (int) { auto tmp = this; operator ++(); return *tmp; }
  ConcatIterator &operator += (inc_t i) { while (i--) { operator ++(); } return *this;  }
  value_type operator * () const { return *it; }
  it_t operator -> () const { return it; }
private:
  const SEQ &s1, &s2;
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
  Concat(const SEQ &s1_,const SEQ &s2_):s1(s1_),s2(s2_) {} 
  const_iterator begin(void) const { return ConcatIterator<SEQ>(s1,s2,s1.begin()); }
  const_iterator end(void) const { return ConcatIterator<SEQ>(s1,s2,s2.end()); }
  size_type size(void) const { return (s1.size() + s2.size()); }
private:
  const SEQ &s1, &s2;
};

void concattest(void) {
  using cv_t = Concat<vector<double>>;
  vector<double> v1 = {1.5,2.5,3.5,4.5};
  vector<double> v2 = {7.0,8.0,9.0,10.0};
  cv_t cv(v1,v2);
  for (typename cv_t::const_iterator it=cv.begin();it!=cv.end();++it)
    cout << *it << ", " << endl;
  cout << "range loop: "; for (auto&& x : cv) cout << x << ", "; cout << endl;
  double s=0.0; for_each(cv.begin(),cv.end(),[&s] (double x) { s += x; });
  cout << "sum = " << s << endl;
  cout << "max = " << *max_element(cv.begin(),cv.end()) << endl;
}

// Template deduction through constructor
template <typename T>
class MyClsTempl {
public:
  using type_t = T;
  MyClsTempl(void) { ptr = nullptr; }
  MyClsTempl(T& x) { ptr = &x; }
  template <typename U>
  T memfn(const T& x,const U& y) const {
    if (x == y) return (x); else return *ptr;
  }
private:
  T *ptr; 
};

// Function template
template <typename ScalarType,typename VectorType>
VectorType saxpy(ScalarType alpha,const VectorType &x,const VectorType &y)
{ return(alpha*x+y); }


// Function type wrappers

double binop(double arg1,double arg2) { return (arg1/arg2); } 

void functionwrapper(void) {
  std::vector<std::function<double(double,double)>> fnvec;
  fnvec.push_back(binop);
  fnvec.push_back([] (double x,double y)  -> double { return y/x; });
  
  for (auto fn : fnvec) { 
    std::cout << fn(3,2) << std::endl;
  }
}

// Initializer list
void initlist(int n = 42) {
  for(long int i: {2,3,4,5,6,7,8,9,n} ) {
    std::cout << i << std::endl;
  }
}

// foreach facility
template <typename T>
void printElements (const T& coll) {
  for (const auto& elem : coll) {
    std::cout << elem << std::endl;
  }
}

// move semantics
class X {
public:
  X(void):n(0),data(nullptr) {}

  // Range initialization
  template <typename It>
  X(It begin,It end):n(0) {
    using value_t = typename It::value_type;
    std::vector<value_t> tmp;
    for(auto it=begin;it!=end;++it,n++) tmp.push_back(*it);
    data = new double [n];
    for(int l=0;l<n;l++) data[l] = tmp[l];
  }

  // Initialization from a collection
  template <typename Coll>
  X(const Coll &c):n(c.size()),data(new double [n]) {
    double *tmp = data;
    for(auto i:c) { *tmp++ = i; }
  }
  // Copy constructor
  X(const X& x):n(x.n),data(new double [n]) {
    cout << "X copy constructor: n = " << n << endl;
    for(int l=0;l<n;l++) data[l] = x.data[l];
  }

  // Move constructor
  X(X &&x):n(x.n),data(x.data) {
    cout << "X move constructor: n = " << n << endl;
    x.n=0; x.data = nullptr; }

  virtual ~X(void) {
    cout << "X destructor: n = " << n << endl;
    if (data != nullptr) delete [] data; }

  X &operator = (const X &x) {
    cout << "X assignment: n = " << x.n << endl;
    if (data != nullptr) delete [] data; 
    n = x.n; data = new double [n];
    for(int l=0;l<n;l++) data[l] = x.data[l];
    return *this;
  }

  X &operator = (X &&x) {
    cout << "X assign & move: n = " << n << ", x.n = " << x.n << endl;
    if (data != nullptr) delete [] data;
    // Stealing the pointers from Rvalue object
    n = x.n; data = x.data;
    x.n = 0; x.data = nullptr;
    return *this;
  }

  friend ostream &operator <<(ostream &o,const X&);
private:
  int n;
  double *data;
};

ostream &operator <<(ostream &o,const X& x) {
  o << '['; for(int l=0;l<x.n;l++) o << x.data[l] << ' ';
  return o << "] ";
}

template<typename T>
std::tuple<T,T,std::vector<T>> extcumsum(const std::vector<T> &v) {
  // Local summation variable captured by reference by lambda function
  T sum{};
  // temporary vector for returning cumulative sum
  std::vector<T> w{};
  // cumulative summation
  std::transform(v.cbegin(),v.cend(),back_inserter(w),
		 [&sum] (T x) { sum += x; return(sum); });
  return(std::tuple<T,T,std::vector<T>>
    (*std::min_element(v.cbegin(),v.cend()),
     *std::max_element(v.cbegin(),v.cend()),
     std::move(w)));
}

void extcumsumdemo(void) {
  // initialize a vector from an initializer list
  std::vector<double> v({1.2,2.3,3.4,4.5,5.6,6.7,7.8});
  // Variables for return values
  double minv,maxv;  // Extremal elements
  std::vector<double> cs; // Cumulative sums
  std::tie(minv,maxv,cs) = extcumsum(v);
  cout << "min = " << minv << ", max = " << maxv << endl;
  cout << "cs = [ "; for(double x: cs) cout << x << ' '; cout << "]" << endl;
}

// Returning several objects
std::tuple<string,X,std::size_t> multireturn(const std::vector<double> &v) {
  std::tuple<string,X,std::size_t> t("tuple",X(v.begin(),v.end()),v.size());
  return(std::move(t));
}

void tupletest_copy(void) {
  cout << "Testing tuples: copy" << endl;
  std::vector<double> v({1.2,2.3,3.4,4.5});
  // Variables for return values
  string s;
  X x{};
  std::size_t sz;

  std::tuple<string,X,std::size_t> t(multireturn(v));
  cout << "X(tuple) = " << std::get<1>(t) << endl;
  std::tie(s,x,sz) = t;
  cout << "s = " << s << ", x = " << x << ", sz = " << sz << endl;
}

void tupletest_move(void) {
  cout << "Testing tuples: move" << endl;
  std::vector<double> v({1.2,2.3,3.4,4.5});
  // Variables for return values
  string s{}; X x{}; std::size_t sz;

  std::tie(s,x,sz) = std::move(multireturn(v));
  cout << "s = " << s << ", x = " << x << ", sz = " << sz << endl;
}

// Demonstration of lambda function and transform algorithm
void lambdademo(void) {
  // initialize a vector from an initializer list
  std::vector<double> v({1.2,2.3,3.4,4.5,5.6,6.7,7.8});
  // A vector of the same length
  std::vector<double> w(v.size());
  // Do cumulative summation of v and store result in w
  double sum = 0;
  std::transform(v.begin(),v.end(),w.begin(),
		 [&sum] (double x) { sum += x; return sum;});
  cout << "sum = " << sum << ", w = [ ";
  for(auto x: w) cout << x << ' '; cout << ']' << endl;
}

int main(int argc,char **argv)
{
  cout << "Testing of C++11 features" << endl;
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " <selection>" << endl;
    return(-1L);
  }
  else {
    const int sel = atoi(argv[1]);
    switch (sel) {
    case 8: { extcumsumdemo(); break; }
    case 7: { lambdademo(); break; }
    case 1: { initlist(); break; }
    case 2: { printElements(std::vector<double>{1.2,2.3,4.5}); break; }
    case 3: {
      std::vector<double> v({1.2,2.3,3.4,4.5});
      X x(v);
      X y(std::move(x));
      cout << "x = " << x << "y = " << y << endl;
      break;
    }
    case 4: { tupletest_copy(); break; } 
    case 5: { tupletest_move(); break; }
    case 6: {
      double x = 3.14;
      float y = 3.14;
      MyClsTempl<int> inst1;
      MyClsTempl<double> inst2(x);
      cout << inst2.memfn(x,y) << endl;
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
    default: { cerr << "Invalid selection" << endl; exit(-1L); }
    }
  }
  return(0);
}
pwd

