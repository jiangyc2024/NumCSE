void storageOrder(int nrows=6,int ncols=7)
{
  cout << "Different matrix storage layouts in Eigen" << endl;
  // Template parameter \texttt{ColMajor} selects column major data layout
  Matrix<double,Dynamic,Dynamic,ColMajor> mcm(nrows,ncols);
  // Template parameter \texttt{RowMajor} selects row major data layout
  Matrix<double,Dynamic,Dynamic,RowMajor> mrm(nrows,ncols);
  // Direct initialization; lazy option: use \texttt{int} as index type
  for (int l=1,i= 0; i< nrows; i++)
    for (int j= 0; j< ncols; j++,l++)
      mcm(i,j) = mrm(i,j) = l;
  
  cout << "Matrix mrm = " << endl << mrm << endl;
  cout << "mcm linear = ";
  for (int l=0;l < mcm.size(); l++) cout << mcm(l) << ',';
  cout << endl;

  cout << "mrm linear = ";
  for (int l=0;l < mrm.size(); l++) cout << mrm(l) << ',';
  cout << endl;
}
