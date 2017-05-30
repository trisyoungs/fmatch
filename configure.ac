#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.58)

# Set program name, version, bug-address and source directory
m4_define([FMATCH3_VERSION],[0.2])
AC_INIT(fmatch3,FMATCH3_VERSION,t.youngs@qub.ac.uk)
AC_CONFIG_SRCDIR([src/main.cpp])

# Check system type and set defines
# Specify a single file in the source tree so we can check that --srcdir points to the right place
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADER(config.h)
AC_CANONICAL_TARGET
MAC_LIBS=""
if test "$target_vendor" = "apple"; then
  AC_DEFINE([IS_MAC], [], [Whether operating system is an Apple])
fi
AC_SUBST(MAC_LIBS)

AM_INIT_AUTOMAKE(fmatch3,FMATCH3_VERSION)

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_PROG_YACC


#
# Additional Configuration Options
#

# Build Directory (used when packaging)
BUILDDIR=""
AC_ARG_WITH([build-dir],[AC_HELP_STRING([--with-build-dir=<dir>],
    [Specify location of build environment when packaging [default=""]])], ,
    [with_build_dir=""])
BUILDDIR="$with_build_dir"

# Install Directory (used when packaging - should be '$PREFIX - $BUILDDIR')
INSTALLDIR="/usr/local"
AC_ARG_WITH([install-dir],[AC_HELP_STRING([--with-install-dir=<dir>],
    [Specify 'clean' instal dir prefix, i.e. without any BUILDDIR [default="/usr/local"]])], ,
    [with_install_dir="/usr/local"])
# If --prefix has been give, use this value by default
if test "$prefix" != "NONE"; then
  INSTALLDIR="$prefix"
fi

if test "$with_install_dir" != "/usr/local"; then
  INSTALLDIR="$with_install_dir"
fi


# Pre-initialise pkg-config  (unless using Qt framework on the Mac)
if test "$with_qt" != "framework"; then
  PKG_PROG_PKG_CONFIG
fi


VERSION=FMATCH3_VERSION
AC_SUBST(BUILDDIR)
AC_SUBST(INSTALLDIR)
AC_SUBST(VERSION)


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_STRUCT_TM


# Checks for library functions.
AC_FUNC_SELECT_ARGTYPES
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor pow select sqrt strchr strrchr strstr])


# Perform Makefile substitutions
AC_LD_FLAGS=" -lreadline"
AC_SUBST(AC_LD_FLAGS)


# Define Makefiles in AC_OUTPUT 
AC_OUTPUT(Makefile src/Makefile)
