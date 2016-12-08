//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#pragma once

// BEGIN: Advanced demonstration code
// The following functions demonstrate an implementation
// of the taylor expansion based computation of sinh.
// The implementation is done efficiently using SFINAE.
// The implementation is stable using a taylor expansion.
// The implementation can be tuned to be as precise as wanted.

/* \brief Recursively compiute factorial.
 * Factorial is computed at compile time to maximize efficiency.
 * \param[in] n an integer $n \geq 0$.
 * \return $n!$
 */
inline constexpr unsigned int factorial(unsigned int n) {
    // If $n > 0$ then $n! = n \cdot (n-1)!$
    // else
    return n > 0 ? n * factorial( n-1 ) : 1;
}

/* This is a dummy class to avoid necessity of forward
 * declare pow(). This is an hack, maybe there is
 * a better way.
*/
class _pow {
public:

    /* \brief Compute the 1-th power of a number (i.e identity)
     * \tparam N SFINAED to accept only 1
     * \param x Number for which compute $x^N$
     * \return $x^1$
     */
    template <unsigned int N>
    inline static double pow(double x,
                      // Advanced C++ template metaprogramming using SFINAE:
                      // Next line is just a dummy line: if N != 1
                      // the scruct std::enable_if<false, bool> has no
                      // type ::type, therefore compilation will fail
                      // when compilation fails, error is not thrown
                      // if another valid template is found
                      typename  std::enable_if<N == 1, bool>::type = false
                      ) {
        // Identity
        return x;
    }

    /* \brief Compute the N-th power of a number, N even
     * \tparam N SFINAED to accept only even numbers
     * \param x Number for which compute $x^N$
     * \return $x^N$
     */
    template <unsigned int N>
    inline static double pow(double x,
                      // Advanced C++ template metaprogramming using SFINAE:
                      // Next line is just a dummy line: if N == 1 or N odd
                      // the scruct std::enable_if<false, bool> has no
                      // type ::type, therefore compilation will fail
                      // when compilation fails, error is not thrown
                      // if another valid template is found
                      typename std::enable_if<N != 0 && (N % 2 == 0), bool>::type
                      = false) {
        // $x^N = x^{N/2}x^{N/2}$
        double y = pow<N/2>(x);
        return y*y;
    }

    /* \brief Compute the N-th power of a number, N odd
     * \tparam N SFINAED to accept only odd numbers
     * \param x Number for which compute $x^N$
     * \return $x^N$
     */
    template <unsigned int N>
    inline static double pow(double x,
                      // Advanced C++ template metaprogramming using SFINAE:
                      // Next line is just a dummy line: if N == 1 or N even
                      // the scruct std::enable_if<false, bool> has no
                      // type ::type, therefore compilation will fail
                      // when compilation fails, error is not thrown
                      // if another valid template is found
                      typename std::enable_if<N != 1 && (N % 2 != 0), bool>::type
                      = false) {
        // $x^N = x^{N/2-1}x^{N/2-1}x$
        double y = pow<(N-1)/2>(x);
        return x*y*y;
    }
};

/* \brief Compute sinh(x) in an efficient way.
 * Use Taylor expansion with length 0.
 * Dummy function to terminate template recursion.
 * \tparam N SFINAED to allow only 0
 * \param x where to evaluate sinh
 * \return sinh(x)
*/
template <unsigned int N,
          // Dummy template parameter exploiting SFINAE
          // if N != 0, type does not exist and
          // template instantiation fails, error is not thrown
          typename std::enable_if<N == 0, bool>::type = true >
inline double taylor_sinh (double x) {
    // Taylor expansion with empty sum is 0
    return 0;
}

/* \brief Compute sinh(x) in an efficient way.
 * Use Taylor expansion with length $N$.
 * \tparam N SFINAED to allow only $N > 0$.
 * \param x Where to evaluate sinh
 * \return sinh(x)
*/
template <unsigned int N,
          // Dummy template parameter exploiting SFINAE
          // if N == 0, type does not exist and
          // template instantiation fails, error is not thrown
          typename std::enable_if<(N > 0), bool>::type = false >
inline double taylor_sinh(double x) {
    // Use taylor expansion with $N-1$ terms to compute
    // the one with $N$ terms
    return taylor_sinh<N-1>(x) + _pow::pow<2*N-1>(x) / factorial(2*N-1);
}

template <unsigned int N>
inline double error_bound(double x) {
    return std::exp(x) * _pow::pow<2*N>(x)  / factorial(2*N+1);
}

// END: Advanced demonstration code
