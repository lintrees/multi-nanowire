#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([multi-nanowire], [0.1], [shaojw@student.sysu.edu.cn])
AC_CONFIG_SRCDIR([src/charge_distribution.cpp])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
#AC_PROG_CC
AC_LANG([C++])

# Checks for libraries.
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
CXXFLAGS="$CXXFLAGS -std=c++0x"
AC_TYPE_SIZE_T
AC_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
AC_C_INLINE
if test "$ac_cv_c_inline" != no ; then
AC_DEFINE([HAVE_INLINE], [1], [Define to 1 if use inline function in  the `gsl' library])
#AC_SUBST([HAVE_INLINE])
fi

# Checks for library functions.

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
