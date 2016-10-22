#include <iostream>
#include <list>
#include <ctime>
#include <random>
#include <cmath>

#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace Eigen;

MatrixXd shape_ident_matrix(const MatrixXd & X) {
        assert(X.rows() == 2 && "X must have 2 rows!");
        unsigned n = X.cols();

        MatrixXd B = MatrixXd::Zero(2*n, 4);

        for(unsigned int row = 0; row < n; ++row) {
            B(2*row,  0) = X(0,row);
            B(2*row,  1) = X(1,row);
            B(2*row+1,2) = X(0,row);
            B(2*row+1,3) = X(1,row);
        }

        return B;
}

double solve_lsq(const MatrixXd & X,
                 const MatrixXd & P,
                 MatrixXd & A) {
    MatrixXd B = shape_ident_matrix(X);

    A = Map<MatrixXd>(
                MatrixXd((B.transpose() * B).ldlt().solve(B.transpose() *
                                                 Map<const MatrixXd>(P.data(), 30, 1)
                                                 )).data(),
                2, 2).transpose();

    return (Map<const MatrixXd>(P.data(), 2, 15) - A*X).norm();

}

enum Shape { Tree, Star };

Shape identify(const MatrixXd Xtree,
               const MatrixXd Xstar,
               const MatrixXd & P,
               MatrixXd & A) {
    MatrixXd Atree, Astar;

    double err_tree = solve_lsq(Xtree, P, Atree);
    double err_star = solve_lsq(Xstar, P, Astar);

    std::cout << "Tree residual norm: " << err_tree << std::endl
              << "Star residual norm: " << err_star << std::endl;

    if(err_star >= err_tree) {
        std::cout << "Points appear to define a tree!" << std::endl;
        A = Atree;
        return Tree;
    }
    else {
        std::cout << "Points appear to define a star!" << std::endl;
        A = Astar;
        return Star;
    }
}


#if INTERNAL
MatrixXd generate_noisy_shape(const MatrixXd & X, double eps) {
    MatrixXd Xnoisy;
    Xnoisy.resizeLike(X);

    MatrixXd A = MatrixXd::Random(2,2);

    for(unsigned int i = 0; i < X.cols(); ++i) {
        Xnoisy.col(i) = A*X.col(i) + eps*VectorXd::Random(2);
    }

    return Xnoisy;

}
#endif // INTERNAL

MatrixXd transform(const MatrixXd & X,
                   const MatrixXd & A) {
    MatrixXd Xret;
    Xret.resizeLike(X);

    for(unsigned int i = 0; i < X.cols(); ++i) {
        Xret.col(i) = A*X.col(i);
    }

    return Xret;
}

void plot(const MatrixXd & Xshape,
          const MatrixXd & X,
          const std::string & title,
          const std::string & name,
          const MatrixXd & A = MatrixXd::Identity(2,2)) {
    MatrixXd Xtransform = transform(X, A);

    mgl::Figure fig;
    fig.xlabel("x");
    fig.ylabel("y");
    fig.plot(Xtransform.row(0), Xtransform.row(1), " *b").label("Noisy points");
    if(Xshape.rows() == 2)
        fig.plot(Xshape.row(0), Xshape.row(1), "-r").label("Training shape");
    fig.title(title);
    fig.legend();
    fig.save(name + ".eps");
}

int main(int argc, char **argv) {
    const unsigned int n = 15;
    MatrixXd Xtree(2, n);
    Xtree << 1, 4,  1,  1,  -1, -1, -4, -1, -3, -1, -2, 0, 2, 1, 3,
             0, -3, -3, -4, -4, -3, -3, 0,  0,  2,  2,  4, 2, 2, 0;

    MatrixXd Xstar(2, n);
    Xstar << // x-coords
             1.00000,  0.40451,  0.30902, -0.15451,  -0.80902,
             -0.50000, -0.80902, -0.15451, 0.30902, 0.40451,
             // x-coords midpoint
             0.60301, 0.18634, -0.48784, -0.48784, 0.18634,
             // y-coords
             0.00000, 0.29389,  0.95106,  0.47553,  0.58779,
             0.00000, -0.58779, -0.47553, -0.95106, -0.29389,
             // x-coords midpoint
             0.00000, 5.7349e-01, 3.5444e-01, -3.5444e-01, -5.7349e-01;

    MatrixXd _Xtree(2,n+1);
    _Xtree << Xtree, Xtree.leftCols<1>();
    MatrixXd _Xstar(2, n-5+1);
    _Xstar << Xstar.leftCols<(n-5)>(), Xstar.leftCols<1>();;

    auto points = [&_Xtree, &_Xstar] (Shape s) {
        return s == Tree ? _Xtree : _Xstar;
    };

#if INTERNAL
    if(argc > 1) {
        srand(time(nullptr));

        MatrixXd X;
        double eps = 10e-3;
        if(argc > 2) {
            eps = std::stof(argv[2]);
        }

        if( std::stoi(argv[1]) == 1)
            X = generate_noisy_shape(Xtree, eps);
        else
            X = generate_noisy_shape(Xstar, eps);

        for(unsigned int i = 0; i < n; ++i) {
            std::cout << X(0, i) << ", ";
        }
        for(unsigned int i = 0; i < n; ++i) {
            std::cout << X(1, i) << ", ";
        }

        return 0;
    }
#endif // INTERNAL

    {
        MatrixXd P(2, n);
        P << -0.370407, 0.737504, 1.85296, 2.61666, 3.36683,
            2.61022, 3.74738, 0.370185, 1.11543, -1.10877,
            -0.745235, -2.97745, -2.23456, -1.86617, -1.12343,
            0.02186, -1.1121, -1.17265, -1.58382, -1.61935,
            -1.20599, -1.26846, -0.0132163, -0.0618339, 0.784893,
            0.765241, 1.60338, 0.825417, 0.808752, 0.0582148;

       std::cout << "****************** Set 1 ******************"
                  << std::endl;
       MatrixXd A;
       Shape s = identify(Xtree, Xstar, P, A);
       plot(transform(points(s), A), P,
            "Set 1: original points",
            "set-1-original");
       plot(points(s), P,
            "Set 1: points with best match",
            "set-1-transformed", A.inverse());
    }
    {
        MatrixXd P(2, n);
        P << 0.00384823, -0.594978, -0.453248, -0.537803, -0.253519,
                -0.137096, 0.0778092, 0.0626776, 0.196791, 0.272765,
                0.330656, 0.356591, 0.104327, 0.0293895, -0.187963,
                0.251845, -0.218366, -0.674227, -0.8617, -1.33682,
                -1.10144, -1.54197, -0.166215, -0.526587, 0.309162,
                0.262935, 1.17761, 0.857399, 0.719816, 0.562366;

       std::cout << "****************** Set 2 ******************"
                  << std::endl;
       MatrixXd A;
       Shape s = identify(Xtree, Xstar, P, A);
       plot(transform(points(s), A), P,
            "Set 2: original points",
            "set-2-original");
       plot(points(s), P,
            "Set 2: points with best match",
            "set-2-transformed", A.inverse());
    }
    {
        MatrixXd P(2, n);
        P << -0.867156, -0.11333, 0.252925, 0.509669, 1.1114,
                0.526225, 0.42909, -0.110557, -0.846383, -0.495474,
                -0.598265, 0.184648, 0.694457, 0.259439, -0.558467,
                -0.65055, -0.499286, -0.751277, -0.201265, 0.292899,
                0.346837, 0.807498, 0.253455, 0.240331, -0.121647,
                -0.380944, -0.477591, 0.144312, 0.554164, 0.135438;

        std::cout << "****************** Set 2 ******************"
                  << std::endl;
        MatrixXd A;
        Shape s = identify(Xtree, Xstar, P, A);
        plot(transform(points(s), A), P,
             "Set 3: original points",
             "set-3-original");
        plot(points(s), P,
             "Set 3: points with best match",
             "set-3-transformed", A.inverse());
    }


    return 0;
}
