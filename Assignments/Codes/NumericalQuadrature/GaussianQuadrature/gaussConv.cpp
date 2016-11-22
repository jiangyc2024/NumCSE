# include <cmath>
# include <vector>
# include <Eigen/Dense>

# include <figure/figure.hpp>

# include "gaussquad.hpp"

/*!
 * \brief integrate Compute integral given quadrature rule.
 * Compute the integral $\int_a^b f(x) dx \approx \sum_{i=1}^n w_i f(c_i)$
 * for given quadrature rule $\{(w_i, x_i)\}_{i=1}^n$
 * \tparam Function A function object with operator()(double)
 * \param qr A QuadRule object passing nodes and weights
 * \param f A function handle passing the integrand
 * \return Value of integral of $f$ using quadrature rule Q
 */
template <class Function>
double integrate(const QuadRule& qr, const Function& f) {
    double I = 0;
    for (unsigned i = 0; i < qr.weights.size(); ++i) {
        I += qr.weights(i) * f(qr.nodes(i));
    }
    return I;
}

/*!
 * \brief gaussConv Approximte the integral $\int_{-1}^1 \arcsin(x) f(x) dx$
 *
 * \tparam Function A function object with operator()(double)
 * \param fh Will pass the integrand
 * \param I_ex Exact value of integral
 * \return Value of integral
 */
/* SAM_LISTING_BEGIN_1 */
template <class Function>
void gaussConv(const Function& fh, const double I_ex) {

    std::vector<double> evals, // used to save no. of quad nodes
            error; // used to save the error
#if SOLUTION
    // Build Integrand
    auto f = [fh](double x) {
        return std::asin(x)*fh(x);
    };

    const unsigned N = 50; // Max. no. of nodes

    for (unsigned n = 1; n <= N; ++n) {
        QuadRule qr = gaussquad(n); // Create quadrature rule

        double I = integrate(qr, f); // Compute integral

        evals.push_back(n); // Save no. of quadrature nodes
        error.push_back(std::abs(I - I_ex)); // Save error
    }
#else // TEMPLATE
    // TODO: compute vector of no. of nodes and of error
#endif

    // Create convergence plots
    mgl::Figure fig;
    fig.title("Gauss quadrature convergence");
    fig.setlog(true, true); // log-log scaling
    fig.plot(evals, error, " +r").label("Error"); // plot error
    fig.fplot("x^(-3)", "k--").label("O(n^{-3})"); // reference line
    fig.xlabel("No. of quadrature nodes");
    fig.ylabel("|Error|");
    fig.legend();
    fig.save("GaussConv");
}
/* SAM_LISTING_END_1 */

/*!
 * \brief gaussConv Approximte the integral $\int_{-1}^1 \arcsin(x) f(x) dx$
 * Ensures that convergenge is expoenential using appropriate transformation,
 * provided $f$ is a smooth function.
 * \tparam Function A function object with operator()(double)
 * \param fh Will pass the integrand
 * \param I_ex Exact value of integral
 * \return Value of integral
 */
/* SAM_LISTING_BEGIN_2 */
template <class Function>
void gaussConvCV(const Function& f, const double I_ex) {
    std::vector<double> evals, // Used to save no. of quad nodes
                        error; // Used to save the error
#if SOLUTION
    // Transform integrand
    auto g = [f] (double x) {
        return x*f(std::sin(x))*std::cos(x);
    };

    const unsigned N = 50; // Max. no. of nodes

    for (unsigned n = 1; n <= N; ++n) {
        QuadRule qr = gaussquad(n); // Obtain quadrature rule

        // Transform nodes and weights to new interval
        Eigen::VectorXd w = qr.weights * M_PI/2;
        Eigen::VectorXd c = ( -M_PI/2 + M_PI/2*(qr.nodes.array() + 1) ).matrix();

        // Evaluate $g$ at quadrature nodes $c$
        Eigen::VectorXd gc = c.unaryExpr(g);

        // Same as $I = \sum_{i=1}^n w_i g(c_i)$
        double I = w.dot(gc);

        evals.push_back(n); // Save no. of quadrature nodes
        error.push_back(std::abs(I - I_ex)); // Save error

    }
#else // TEMPLATE
    // TODO: compute vector of no. of nodes and of error (exponential convergence)
#endif

    // create convergence plots
    mgl::Figure fig;
    fig.title("Gauss quadrature convergence");
    fig.setlog(false, true); // lin-log scaling
    fig.plot(evals, error, " +b").label("Error"); // plot error
    fig.xlabel("No. of quadrature nodes");
    fig.ylabel("|Error|");
    fig.legend();
    fig.save("GaussConvCV");
}
/* SAM_LISTING_END_2 */

int main() {
    // "Exact" value of integral
    const double I_ex = 0.870267525725852642;

    // $f(x) = \sinh x$
    std::function<double(double)> fh = [](double x) {
        return std::sinh(x);
    };

    // PART 1
    gaussConv(fh, I_ex);

    // PART 2
    gaussConvCV(f, I_ex);

    return 0;
}
