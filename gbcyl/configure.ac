#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(gbcyl, 6.4.0, jouni.karjalainen@oulu.fi)
AM_INIT_AUTOMAKE([-Wall -Werror foreign tar-pax])
AC_PROG_FC
AC_PROG_CC
# Checks for programs.
AC_PROG_RANLIB
AC_PROG_LIBTOOL

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AX_F90_MODULE_EXTENSION
AC_SUBST(MODEXT, [$ax_f90_modext])

# Checks for library functions.

# Doxygen support:
DX_HTML_FEATURE(ON)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN([$PACKAGE_NAME], [doxygen.cfg], [doxygen-doc])


AC_CONFIG_FILES([
Makefile
mt_stream_f90-1.11/Makefile
ftnunit/Makefile
src/Makefile
data/Makefile
])
AC_OUTPUT

