#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(ptgbcyl, 1.0.0, jouni.m.karjalainen@gmail.com)
AM_INIT_AUTOMAKE([-Wall -Werror subdir-objects foreign tar-pax])
AC_PROG_FC
AC_PROG_CC
AC_PROG_F77
# Checks for programs.
AC_CONFIG_MACRO_DIR([m4])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AX_F90_MODULE_EXTENSION
AC_SUBST(MODEXT, [$ax_cv_f90_modext])

# Checks for library functions.

# Doxygen support:
DX_HTML_FEATURE(ON)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN([$PACKAGE_NAME], doxygen.cfg, doxygen-doc)

# Python
AM_PATH_PYTHON([2.6],,[:])

AC_ARG_ENABLE([pfunit],
[  --enable-pfunit Turn on tests with pFUnit],
[case "${enableval}" in
  yes) pfunit=true ;;
  no)  pfunit=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-pfunit]) ;;
esac], [pfunit=false])
AM_CONDITIONAL([PFUNIT], [test x$pfunit = xtrue])

AC_ARG_ENABLE([mttests],
[  --enable-mttests Turn on mt_stream_f90 tests.],
[case "${enableval}" in
  yes) mttests=true ;;
  no)  mttests=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-mttests]) ;;
esac], [mttests=false])
AM_CONDITIONAL([MT_TESTS], [test x$mttests = xtrue])

AC_ARG_ENABLE([claftests],
[ --enable-claftests Turn on claf90 tests.],
[case "${enableval}" in
  yes) claftests=true ;;
  no) claftests=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-claftests]) ;;
esac], [claftests=false])
AM_CONDITIONAL([CLAF_TESTS], [test x$claftests = xtrue])

AC_CONFIG_FILES([
Makefile
claf90/Makefile
mt_stream_f90-1.11/Makefile
cephes/Makefile
BLAS/Makefile
lapack/lapack_routine/Makefile
lapack/Makefile
src/Makefile
Examples/Makefile
Examples/LJ_and_GB_in_R7/Makefile
tests/Makefile
])
AC_OUTPUT

