#ifndef ARROWMATRIX_FUNCTIONS_HPP
#define ARROWMATRIX_FUNCTIONS_HPP

#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <vector>

#include "timer.h"
#include "matplotlibcpp.h"
#include "plot.hpp"

using namespace Eigen;
namespace plt = matplotlibcpp;


/* @brief Build an "arrow matrix" and compute A*A*y
 * Given ve[conf.yml](/cx_project_file/sXBQr9N7BqZNQ356r?download=true)ctors $a$ and $d$, returns A*A*x in $y$, where A is built from a, d
 * @param[in] d An n-dimensional vector
 * @param[in] a An n-dimensional vector
 * @param[in] x An n-dimensional vector
 * @param[out] y The vector y = A*A*x
 */
/* SAM_LISTING_BEGIN_0 */
void arrow_matrix_2_times_x(const VectorXd &d, const VectorXd &a,
                            const VectorXd &x, VectorXd &y) {
  assert(d.size() == a.size() && a.size() == x.size() &&
         "Vector size must be the same!");
  int n = d.size();

  // In this lines, we extract the blocks used to construct the matrix A.
  VectorXd d_head = d.head(n-1);
  VectorXd a_head = a.head(n-1);
  MatrixXd d_diag = d_head.asDiagonal();

  MatrixXd A(n,n);

  // We build the matrix A using the "comma initialization": each expression separated
  // by a comma is a "block" of the matrix we are building.
  // d\_diag is the top left (n-1)x(n-1) block
  // a\_head is the top right vertical vector
  // a\_head.transpose() is the bottom left horizontal vector
  // d(n-1) is a single element (a 1x1 matrix), on the bottom right corner
  // This is how the matrix looks like:
  // A = | D   | a      |
  //     |-----+--------|
  //     | a\^T | d(n-1) |
  A << d_diag,             a_head,
       a_head.transpose(), d(n-1);

  y = A*A*x;
}
/* SAM_LISTING_END_0 */


/* @brief Build an "arrow matrix"
 * Given vectors $a$ and $b$, returns A*A*x in $y$, where A is build from a,d
 * @param[in] d An n-dimensional vector
 * @param[in] a An n-dimensional vector
 * @param[in] x An n-dimensional vector
 * @param[out] y The vector y = A*A*x
 */
/* SAM_LISTING_BEGIN_1 */
void efficient_arrow_matrix_2_times_x(const VectorXd &d,
                                      const VectorXd &a,
                                      const VectorXd &x,
                                      VectorXd &y) {
  assert(d.size() == a.size() && a.size() == x.size() &&
       "Vector size must be the same!");
  int n = d.size();

  // TO DO: (2-1.c) Implement an efficient version of arrow_matrix_2_times_x.
  // Hint: Notice that performing two matrix vector multiplications is less
  // expensive than a matrix-matrix multiplication.
  // Therefore, as first step, we need a way to efficiently compute A*x
  
  // START
  
  // END
}
/* SAM_LISTING_END_1 */


/* \brief Compute the runtime of arrow matrix multiplication.
 * Repeat tests 10 times, and output the minimal runtime
 * amongst all times. Test both the inefficient and the efficient
 * versions.
*/
/* SAM_LISTING_BEGIN_3 */
void runtime_arrow_matrix() {
  //Memory allocation for plot
  std::vector<double> vec_size;
  std::vector<double> elap_time, elap_time_eff;


  //header for result print out 
  std::cout << std::setw(8) << "n"
            << std::setw(15) << "original"
            << std::setw(15) << "efficient"
            << std::endl;


            
            
  for (unsigned n = 8; n <= 128; n *= 2) {

    //save vector size (for plot)
    vec_size.push_back(n);

    // Number of repetitions
    unsigned int repeats = 10;

    Timer timer, timer_eff;
    // Repeat test many times
    for (unsigned int r = 0; r < repeats; ++r) {
      // Create test input using random vectors
      Eigen::VectorXd a = Eigen::VectorXd::Random(n),
                      d = Eigen::VectorXd::Random(n),
                      x = Eigen::VectorXd::Random(n),
                      y;

      // Compute times for original implementation
      timer.start();
      arrow_matrix_2_times_x(d, a, x, y);
      timer.stop();

       
      // TO DO: (2-1.e) Compute times for efficient implementation
      // START
      
      // END

    }



    // Print results (for grading): inefficient 
    std::cout << std::setw(8) << n
              << std::scientific << std::setprecision(3)
              << std::setw(15) << timer.min()
              << std::setw(15) << timer_eff.min()
              << std::endl;        


    //time needed
    elap_time.push_back(timer.min());
    elap_time_eff.push_back(timer_eff.min());
  }

  /* DO NOT CHANGE */ 
  //create plot
  plot(vec_size, elap_time, elap_time_eff, "./cx_out/text.png");
}
/* SAM_LISTING_END_3 */


#endif
