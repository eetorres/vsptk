#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([vsptk], [1.0], [eetorres@gmail.com])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/vsptk.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor sqrt strrchr strstr])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
