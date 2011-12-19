AC_INIT([TreeMix], [0.1], [joepickrell@gmail.com], [treemix], [http://treemix.googlecode.com/])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_CONFIG_HEADERS([config.h])
AC_CHECK_LIB([z],[gzopen],,AC_MSG_ERROR(could not find zlib))
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,AC_MSG_ERROR(could not find BLAS))
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,AC_MSG_ERROR(could not find GSL))
AC_PROG_CXX
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
