AC_PREREQ([2.69])
AC_INIT([sage],[6.0],[http://groups.google.com/group/sage-support])
AC_CONFIG_FILES([Makefile])

m4_include([m4/ax_compare_version.m4])

gcc_message="You can try to install a recent version of GCC and rerun this script. If you want to try to build sage with the compiler that is available, then call this script with the option --disable-gcc-check (not recommended). If you want us to build a recent GCC which it then uses during the build process, then call this script with the option --with-build-gcc."

AC_PROG_CC
AS_IF([test x"$CC" = x],[AC_MSG_ERROR([A C compiler is required to build sage but none was found. Please install a C compiler.])])

AC_ARG_ENABLE([gcc-check], AS_HELP_STRING([--disable-gcc-check], [Disable gcc version checks]))
AC_ARG_WITH([build-gcc], AS_HELP_STRING([--with-build-gcc], [Build a recent version of GCC as soon as portage is usable (this option implies --disable-gcc-check)]))
AS_IF([test x"$with_build_gcc" = xyes], [enable_gcc_check=no])
AS_IF([test x"$enable_gcc_check" != xno], [
	AS_IF([test x"$GCC" = xyes],,[AC_MSG_ERROR([The GNU C compiler is recommended to build sage but was not found. $gcc_message])])

	gcc_version="`$CC -dumpversion`"
	AX_COMPARE_VERSION([$gcc_version],[lt],[4.4],[AC_MSG_ERROR([Your version ($gcc_version) of the GNU C compiler is too old to build sage. $gcc_message])])
	AX_COMPARE_VERSION([$gcc_version],[eq0],[4.6.0],[AC_MSG_ERROR([Your version ($gcc_version) of the GNU C compiler is known not to work when building sage. $gcc_message])])
	AX_COMPARE_VERSION([$gcc_version],[eq0],[4.6.1],[AC_MSG_ERROR([Your version ($gcc_version) of the GNU C compiler is known not to work when building sage. $gcc_message])])
])

AC_ARG_WITH([python], AS_HELP_STRING([--with-python], [Build Python to be able to run portage]))
AM_PATH_PYTHON([2.7],,[
	AS_IF([test x"$with_python" = xyes], [
		AC_MSG_ERROR([The sage build system uses portage which requires the installed Python to be at least version 2.7. Your system does not have a sufficiently recent version of Python installed. You can try to install a recent version of Python and rerun this script. If you want us to build a recent version of python to run portage, call this script with the option --with-python])
	])
])

AC_OUTPUT()
