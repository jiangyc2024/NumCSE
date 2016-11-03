///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Till Ehrengruber <tille@ehtz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

// return name of eigen's fft backend
inline std::string eigen_fft_backend_name() {
#ifdef EIGEN_FFTW_DEFAULT
	// FFTW: faster, GPL -- incompatible with Eigen in LGPL form, bigger code size
	return "FFTW";
#elif defined EIGEN_MKL_DEFAULT
	// TODO 
	// intel Math Kernel Library: fastest, commercial -- may be incompatible with Eigen in GPL form
	return "Intel Math Kernel Library";
#else
	// ei_kissfft_impl:  small, free, reasonably efficient default, derived from kissfft
	return "Kiss FFT";
#endif
}
