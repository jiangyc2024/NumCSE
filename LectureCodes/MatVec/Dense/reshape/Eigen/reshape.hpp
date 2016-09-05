template<typename MatType>
void reshapetest(MatType &M)
{
  using index_t = typename MatType::Index;
  using entry_t = typename MatType::Scalar;
  const index_t nsize(M.size());

  // reshaping possible only for matrices with an even number of entries
  if ((nsize %2) == 0) {
    entry_t *Mdat = M.data();
    Map<Eigen::Matrix<entry_t,Dynamic,Dynamic>> R(Mdat,2,nsize/2);
    cout << "Matrix M = " << endl << M << endl;
    cout << "reshaped to " << R.rows() << 'x' << R.cols()
	 <<  " = " << endl << R << endl;
    // Modifying R \textbf{affects M}, because they share the data space !
    R *= -1.5;
    cout << "Matrix M = " << endl << M << endl;
  }
}