///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Till Ehrengruber <tille@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
///
/// This file contains compatibility fixes for gcc versions below 5 which
///  have limited c++11 supported.
//////////////////////////////////////////////////////////////////////////
#include <iostream>

#ifdef __GNUC__ // bug in gcc < 5 has incomplete iostream
// https://gcc.gnu.org/viewcvs/gcc?view=revision&revision=215952
#ifndef __clang__ // clang pretends to be gcc, but is working fine
#if __cplusplus >= 201103L && __GNUG__ < 5
namespace std {

inline ios_base&
defaultfloat(ios_base& __base)
{
    __base.unsetf(ios_base::floatfield);
    return __base;
}

inline ios_base&
hexfloat(ios_base& __base)
{
    __base.setf(ios_base::fixed | ios_base::scientific,
                ios_base::floatfield);
    return __base;
}

}
#endif
#endif
#endif
