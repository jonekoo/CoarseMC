#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(gbcyl, 1.1, jounik@mail.student.oulu.fi)
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_FC
# Checks for programs.

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AX_F90_MODULE_EXTENSION
AC_SUBST(MODEXT, [$ax_f90_modext])

# Checks for library functions.

AC_CONFIG_FILES([
Makefile
src/Makefile
data/Makefile
])
AC_OUTPUT
