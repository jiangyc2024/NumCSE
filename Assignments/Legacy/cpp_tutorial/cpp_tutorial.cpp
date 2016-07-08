#include <Eigen/Dense>

#include <iostream>
#include <vector>
#include <functional>

template <typename T> struct hello {
    T operator*() { return "Hello"; }
    static const char e = '!';
    template <char c> char put_char() { return c; }
};

int incr(int a) { return a+1; }

int Id(int a) { std::cout << "int\n"; return a; }
int Id(double a) { std::cout << "double\n"; return a; }
int Id() { std::cout << "void\n"; return 8; }
// error: ambiguating new declaration:
// double Id() { return 8; } 
// double Id(int a = 9) { return 8; }

typedef std::complex<double> Complex;
bool operator<(Complex a, Complex b) {
    return a.real() < b.real();
}

namespace MySpace {
    int i = 9;
}


//// OOP: class and struct
struct myStruct {
    int     some_int;
    double  some_double;
    void    incr() { some_int++; } // function
};
class myClass {
public:
    //// Constructors
    myClass() { I = 9; }; // default constructor
    myClass(int I_) : I(I_) {} // any constructor
    myClass(const myClass & other) // copy constructor
    : I(other.I) {} // special syntax I = other.I
    int        I;
    myStruct   S; // contains a struct
    
    //// Const member
    int nonconst_function() {
        // something
        I = 8; // fine
    }
    int const_function() const {
        // something
        //             I = 8; // is illegal: modifies myClass
    }
    
    //// Static member
    static int static_I;
    void incr() { static_I++; }
    
    //// Static function
    static int function() {
        int a = static_I; // valid
//         int b = I; // is illegal: cannot see non static
        return 8;
    }
    
    //// Operator
    int operator*(const myClass & other) const {
        return 5*I + other.I;
    }
    int weird_mult(const myClass & other) const {
        return 5*I + other.I;
    }
};

int myClass::static_I = 9;

//// OOP: class namespaces
struct N {
    class N_c { // nested class
        
    };
    typedef double N_t;
    static int N_i;
    void N_f(); // function declaration
};

int N::N_i;
void N::N_f() { ++N_i; } // definition of function N_f

//// Tempate
double multAB(double A, double B) { return A*B; }

template <class whatever>
whatever multAB(const whatever & A,
                const whatever & B) { return A*B; }


template <class T1, class T2> // s.t. can mix T1/T2
T1 min(const T1 & a, const T2 & b) {
    return a < b ? a : b; // ternary operator
}

//// Class template
template <class T>
struct myComplex {
    T Re, Im; // real and imaginary part
    
    myComplex(T Re_, T Im_) : Re(Re_), Im(Im_) {};
};

template <class T1, class T2 = int> // T2 default int
class myComplex2 {
    T1 Re; // real part has type T1
    T2 Im; // imaginary part part has type T2
public:
    template <class U>
    U sum() { return Re + Im; }
};

//// Typename
template <class T1, class T2>
void templated_function() {
    int x;
    T1::t * x; // multiplication t * x
    typename T2::t * y; // declaration of y
}

struct t_var {
    static int t; // t is a variable
};

struct t_class {
    class t { // t is a class name
        // stuff...
    };
};
//// template keyword 2nd meaning

template <typename T>
struct tStruct {
    template <typename U>
    U tMemberFunc(void) { return U(); }
};

template <typename T>
void tFunc() {
    T s;
//     int a = s.myMemberFunc<int>(); // error
    int a = s.template tMemberFunc<int>(); // ok
}

/// decltype
template <typename T>
struct declClass {
    T val;
};


/// auto
std::vector< Eigen::MatrixXd > some_function(int a) { return std::vector< Eigen::MatrixXd >(); };
    
int main(void) {
    //// Hello
    std::cout << "Hello world!" << std::endl;

    std::vector<char>  world = {'w','o','r','l','d'};
    hello<std::string> h;
    std::cout << *h << h.put_char<' '>();
    for(auto & word : world) { std::cout << word; }
    std::cout << hello<std::string>::e << "\n";
    
    //// lvalue vs rvalue
    {
        int a;
        a = 9; // a is lvalue, 9 is rvalue
//         9 = a; // invalid!
//         incr(a) = 3; // invalid, f(a) is rvalue
        a = incr(a); // valid
    }
    
    //// Function overload
    Id(1); // prints "int"
    Id(1.); // prints "double"
    Id(); // prints "void"
    
    //// Operator overload
    Complex(3,4) < Complex(3,4); // now makes sense
    
    //// Pointers
    {
        int * a; // ptr. to int, undefined memory location
        int * b = new int; // ptr. to allocated int
        a = b; // a points to the same as b
        *a = 9; // value at the memory pointed by a is 9
        std::cout << *b << "\n"; // b pointed to the same as a
        int c;
        a = &c; // a points to the address of c
        delete b;
    }
    
    //// Casting
    {
        double i = 9; // <=> double i = (double) 9.
    }
    
    //// Namespaces
//     std::cout << i; // error: 'i' was not declared in this scope
    std::cout << MySpace::i << "\n"; // ok
    using namespace MySpace;
    std::cout << i << "\n"; // now it works ok
    int i = 7;
    std::cout << i << "\n"; // i = 7 here
    
    //// Typedef
    typedef std::vector<double> custom_vector;
    custom_vector cv; // cv is std::vector<double>
    
    //// OOP class and scructs
    myClass m; // default construct
    myClass n(1); // construct with ints
    myClass o(n); // copy m to o
    myClass p(); // NOT what you think
    int a(); // function declaration like this
    
    //// this keyword
    struct T {
        int i = 9; // syntax sugar to initialize i = 9
        void weird_recursion(T & myself) {
            if(i > 0) {
                myself.i--;
                weird_recursion(*this);
            }
        }
    };
    T x;
    std::cout << x.i << std::endl; // 9
    x.weird_recursion(x);
    std::cout << x.i << std::endl; // 0
    
    //// Class namespaces
    N::N_c n_c; // <=> n_c is a class N_c
    N::N_t n_d; // <=> double N_d
    N::N_i = 9; // variable of type int
    
    //// OOP: const member
    const myClass mc;
//     mc.nonconst_function(); // illegal
    mc.const_function(); // legal
    
    //// OOP: static variable
    std::cout << m.static_I << " " << m.static_I << "\n";
    m.incr();
    std::cout << m.static_I << " " << m.static_I << "\n";
    
    //// OOP: static function
    int y = myClass::function();
    
    //// OOP: operators
    int ret = m * n; // a custom uncommon operator
    int ret2 = m.weird_mult(n); // is exactly the same
    
    //// OOP: privacy
    class C {
        // Default: private, accessed only within class
        int a; // C.a invalid outside
    public: // can be accessed only by derived class
        int b; // C.b allowed only in child
    protected: // can be accessed by anyone
        int c; // C.c always valid
    } pr; 
    
//     pr.a; // illegal: private
    pr.b;
//     pr.c; // illegal: protected
    
    //// Function templates
    {
        int a = min(8.,9.); // T1 = T2 = double
        double b = min(8,9); // T1 = T2 = int
        double c = min(8,9.); // T1 = int, T2 = double
        // a = b = c = 8
        
        struct S {
            int i;
            S(int i_) : i(i_) { }; // constructor
            operator int() const { return i; } // convert S -> int
            bool operator<(const S & other) const {
                return i > other.i; // swap < for >
            }
        };
        
        // Uses operator< for S, which is actually >
        double d = min(S(8),S(9)); // d = 9
    }
    
    //// Class templates
    myComplex<int>  cmplx1(8,9); // cmplx = 8*i+9
    myComplex<double>  cmplx2(2.3,3); // cmplx = 2.3*i+3
    
    myComplex2<double> cmplx3; // Re = double, Im = int
    myComplex2<double, double> cmplx4; // Re = Im = double
    
    cmplx3.sum<int>(); // must specify U above
    double s1 = cmplx3.sum<double>(); // must specify U above
    
    //// Template typename
    templated_function<t_var, t_class>();
//     templated_function<t_class, t_var>(); // invalid
    
    //// Template keyword
    tFunc<tStruct<int>>();
    
    //// Stl vector
    std::vector<int> vec; // int vec., similar to int *
    std::vector< std::complex<double> > cmplx_vec;
    
    
    std::vector<int> v;
    // ... fill v
    for(int i = 0; i < v.size(); ++i) {
        v.at(i) = i+1;
        // loop over all v
    }
    
    for(std::vector<int>::iterator it = v.begin(); it  != v.end(); ++it) {
        *it = i+1;
        // loop over all v, *it is like v.at(i)
    }
    // range based alternative:
    for(int & c: v) {
        c = i+1;
        // loop over all v, c is like v.at(i)
    }
    
    //// decltype
    declClass<int> m_i;
    declClass<double> m_d;
    decltype(m_i.val) i_type; // <=> int i_type;
    decltype(m_d.val) d_type; // <=> double d_type;
    
    //// auto
    std::vector< Eigen::MatrixXd > vm = some_function(8);
    auto va = some_function(8);
    for(std::vector< Eigen::MatrixXd >::iterator it = va.begin(); it < va.end(); ++it) {
        // some code
    }
    for(auto it = va.begin(); it < va.end(); ++it) {
        // some code
    }
    
    //// Lambda
    auto f = [] (int i) { return i++; }; // using auto
    std::function<int (int)> g = [] (int i) { return i++; }; // manual type
    
    {
        int j = 9, i = 8;
    //     [] () { std::cout << j; }; // illegal!
        // legal, j captured by reference:
        auto f = [&] () { std::cout << j << "\n"; };
        f();
    //     [i,&j] () { i++; j++ }; // illegal for i, legal for j
    }
    
    //// eigen aliasing
    Eigen::VectorXd A(4);
    A << 1,2,3,4;
    A.tail(2) = A.head(2);
//     A[3] = A[2]; A[1] = A[0]; A[2] = A[1];
    std::cout << A << std::endl;
    
//     A = A.transpose(); // aborts during execution
    
    //// Eigen loops
    {
        int n = 9;
        Eigen::MatrixXd A(n,n); // A is Col. Major format
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                A(i,j) = i*j;
            }
        }
        for(int j = 0; j < n; ++j) {
            for(int i = 0; i < n; ++i) {
                A(i,j) = i*j;
            }
        }
    }
}
