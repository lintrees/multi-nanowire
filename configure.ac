#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.68])
AC_INIT([multi-nanowire], [1.0], [shaojw@student.sysu.edu.cn])
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
AC_CHECK_LIB([gmp], [__gmpq_init])
AC_CHECK_LIB([CGAL], [main])



# Checks for header files.
AC_HEADER_ASSERT

# Checks for typedefs, structures, and compiler characteristics.
CXXFLAGS="$CXXFLAGS -pipe -std=c++0x -frounding-math"
AC_TYPE_SIZE_T
AC_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
AC_C_INLINE
if test "$ac_cv_c_inline" != no ; then
AC_DEFINE([HAVE_INLINE], [1], [Define to 1 if use inline function in  the `gsl' library])
fi

AC_DEFINE([SINGLE_SPHERE_CHARGE_APPROX], [], [using single sphere charge approx. when calcing current])

AC_COMPILE_IFELSE(
    [AC_LANG_SOURCE([[constexpr int i(0);]])],
    [],
    [AC_DEFINE([constexpr], [const], [no support for constexpr])]
)

# Checks for library functions.

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
