#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
AC_PREREQ([2.69])
AC_INIT([sage-lib], [none], [sage-devel@googlegroups.com])
AM_INIT_AUTOMAKE([-Wall -Werror -Wno-portability foreign])
AC_CONFIG_SRCDIR([c_lib/configure.ac])

AC_CONFIG_SUBDIRS([bin c_lib sage ext doc])

distdir_subs='distdir: DIST_SUBDIRS+=$(CORE_MODULES)'
AC_SUBST([distdir_subs])
AM_SUBST_NOTMAKE([distdir_subs])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
