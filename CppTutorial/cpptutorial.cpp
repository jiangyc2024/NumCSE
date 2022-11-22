// How to compile this file:
// on Linux (similar on Mac OS X)
// Step 0. Open terminal emulator
// Step 1. Move to directory where source is located
// Step 2. Invoke compiler
//            g++ -Wall -pedantic -std=c++11 cpptutorial.cpp
//         can also use clang++ instead of g++.
// Step 3. Run with ./a.out
// Optional:   if CMake is available and configured, compile using
//         cmake .
//         make
// CMake makes compilations easy and portable, is also used by many IDEs

// The following includes the entire input/output standard Library
#include <iostream>

// This would import all std::XXX functions/objects, so that you do not have to
// use them with std::XXX but simply using XXX instead.
// WARNING: This is fine, but potentially dangerous:
// if you or somebody else defines an ojbect/function with name
// XXX, there may be a conflict of naming.
// The standard library contains common names such as vector, list, min, etc.
// So it is not so uncommon to have conflicting names
// Pro tip: do not use this, you will thank me in the future
//
// using namespace std;
//
// Instead, only import some identifiers from std that are frequently used
using std::cout;
using std::endl;

// Later we will
// using namespace Eigen;

//////////////////////////////
///// CLASSES and STRUCTS
//////////////////////////////

// Similarly as the case of functions, if we want to use a class
// before it is being defined, we need to *declare* the class:
class MyClass; // I Promise to *define* MyClass later

// This is the definition of a struct
struct MyStruct {
// public: // This is implicit in structs
    int a; // All those variables can be acessed from "outside"
    double b;
    bool f;
private: // Can only access from within a struct
    int private_int;
    double private_double;

    // A pointer to a MyClass type
    MyClass* C; // This is why we needed to *declare* myclass
    // This cannot be
    // MyClass C; // This will result in a class of infinite dimension
};

// This is the definition of a class
// Remember: class <=> struct are completely equivalent,
// except in struct members are public by default
// in a class members are private by default
class MyClass {
private: // Default, but doesn't harm to write it anyways
    int a;
    float b;
    char u;
public:
    int public_int; // This can be accessed everywhere

    MyStruct C;
};

//////////////////////////////
// Some function
//////////////////////////////

// This will be a function taking no argument
// and returning nothing
// I promise to define my_function somewhere else
// in the future (either below, or on another file)
void my_function();

// The next function will take an int and
// return the int multiplied by 6
// The integer b is passed by value
// a temporary copy of b is created,
// any modification to b doesn't affect the
// original variable
// IN C++: All objects are passed by value
int multiplyBy6(int b) {
    return 6*b;
}

// The next function  will take a *double*
// and multiply it by 60, also we print "Double"
// to proof that we called this function
// and not the above
double multiplyBy6(double b) {
    cout << "Double" << endl;
    // Multiply and return
    return 60*b;
}

// Take an int, create a "local copy" and increment
// the local copy by 5. This copy is then discarded
// This function does nothing
void increment(int a) {
    a += 5;
}

// Take an int, but do not copy
// Instead create a reference
// b is incremented after this call
void increment_really(int & a) {
    a += 5;
}

// Take an int, and do not copy
// However, promise not to modify a
// Returns 7 times a
int multiplyBy7(const int & a) {
    return 7*a;
}

//////////////////////////////
// main() function is where everything starts
//////////////////////////////

// This is my main function
int main() {
    //////////////////////////////
    // Variables
    //////////////////////////////

    int a = 4;
    double b; // No need to initialize
    // Dangerous, but compiler generally warns you
    bool c;
    char d;

    // Assign value 4
    b = 4;

    //////////////////////////////
    // Output
    //////////////////////////////

    // << is an *operator*, it will "push" string and other types
    // to the standard output
    cout << "This is my string!" << "\n"
              << "This is another string!" << endl // same as "\n"
              << a << " " << b; // can pass int, float, double, etc.
    // Explanation of "<<":
    //    an operator is a special kind of function
    //    for instance the operator * is nothing but a function
    //    instead of writing product(a,b), you write a*b
    //    with "<<" the situation is analogous
    //    we could imagine a function:
    //        std::ostream & print(std::ostream & o, std::string s) {
    //            // write "s" to the stream "o" (somehow)
    //            return o;
    //        }
    //    where ostream represents the "standard output"
    //    e.g. your terminal
    //    but then, writing:
    //        print(print(std::cout, "String 1"), "String 2")
    //    would be very inconvenient, this can be avoided with an operator:
    //        std::cout << "String 1" << "String 2";
    // Explanation of "::":
    //    Tell compiler to find the object "cout" inside the namespace
    //    called "std". A namespace is like a box, with a label, where you
    //    store all of your functions/objects by giving a "prefix" name
    // Explanation of "endl":
    //    "endl" is (almost) another name for "/n", is a bit "safer" (details omitted)

    // The standard lib has many functions with common namespace
    // std::min(5,6);
    // If you "import the std namespace", using namespace std;
    // you can also write:
    // min(5,6); // Without std

    // Call empty function that does nothing
    my_function();

    // Test functions
    cout << "Multiply by 6: " << multiplyBy6(8) << endl;
    // The following is actually calling the double version:
    cout << "Multiply by 6: " << multiplyBy6(8.312) << endl;

    //////////////////////////////
    // References
    //////////////////////////////
    int I = 10;
    increment(I); // here I is copied to a new variable I_prime
    // I has not changed
    cout << "I = " << I << endl; // prints 10
    increment_really(I);
    cout << "I = " << I << endl; // prints 15

    const int J = 9; // I promise that I won't modify J
    multiplyBy7(J); // Valid, passed by reference but primising "const"
    // increment_really(J); // Not valid, will break promise
    increment(J); // Fine: a copy is created

    // Reference type
    int & alias_to_I = I; // I is now another name for alias_to_I
    alias_to_I = 90; // this also changes I
    cout << "Alias_to_I: " << alias_to_I << endl; // prints 90
    cout << "I = " << I << endl; // prints 90

    //////////////////////////////
    //// Pointers
    //////////////////////////////

    int var = 15;
    int * ptr = &var; // &var is the address of var
    // ptr now "points" to var, i.e. it stores the memory location
    // where var is stored

    *ptr = 60; // *prt is "the value pointed by ptr"
    // this is the same as: var = 60;

    cout << "var =  " << var << endl;
    cout << "*ptr = " << *ptr << endl; //< prints the address of var
    cout << "ptr =  " << ptr << endl;

    //////////////////////////////
    // Structures and classes (object oriented programming)
    //////////////////////////////

    // Create a structure, call default constructor
    MyStruct S {};
    S.a;
    // S.private_double; // Cannot access private member here
}

// I promise to define function_b somewhere else
void function_b();

// A calls function_b, and I promised to define function_b somewhere
// this is fine
void function_a() { //NOLINT(misc-no-recursion)
    function_b();
}

// Function_a is defined, no issue here
void function_b() { //NOLINT(misc-no-recursion)
    function_a();
}

// This is my_function definition
// Empty function
void my_function(  ) {
    // Insert something here
}
