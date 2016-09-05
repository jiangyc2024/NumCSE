//! \brief Multiplication of Kronecker product with vector $y = (A \otimes B)x$. Elegant way using reshape
//! WARNING: using Matrix::Map we assume the matrix is in ColMajor format, *beware* you may incur in bugs if matrix is in RowMajor isntead
//! \param[in] A Matrix $m \times n$
//! \param[in] B Matrix $l \times k$
//! \param[in] x Vector of dim $nk$
//! \param[out] y Vector y = kron(A,B)*x of dim $ml$
template <class Matrix, class Vector>
void kronmultv(const Matrix & A, const Matrix & B, const Vector & x, Vector & y)
{
    unsigned int m = A.rows(); unsigned int n = A.cols();
    unsigned int l = B.rows(); unsigned int k = B.cols();
    // 1st matrix mult. computes the products \Blue{$\VB\Vx^{j}$}
    // 2nd matrix mult. compines them linearly with the coefficients of \Blue{$\VA$}
    Matrix t = B * Matrix::Map(x.data(),k,n) * A.transpose();
    y = Matrix::Map(t.data(), m*l, 1);
}