# This looks for the BZLIB libraries. If we find them, we set the Makefile
# conditional HAVE_BZLIB. We as set BZLIB_CPPFLAGS, and BZLIB_LDFLAGS.
# 
# To allow users to build there own copy of BZLIB, we also define
# BUILD_BZLIB

AC_DEFUN([AC_BZLIB],
[
have_bzlib="no"
build_bzlib="no"
AC_ARG_WITH([bzlib],
        AS_HELP_STRING([--with-bzlib@<:@=DIR@:>@], [use bzlib (default is yes if found) - it is possible to specify the root directory for bzlib (optional). You can also specify "build" if you want to build your own local copy. See also THIRDPARTY variable described below. This is only needed if you are building a local copy of boost.]),
        [
    if test "$withval" = "no"; then
        want_bzlib="no"
    elif test "$withval" = "yes"; then
        want_bzlib="yes"
        ac_bzlib_path=""
    elif test "$withval" = "build"; then
        want_bzlib="yes"
        build_bzlib="yes"
        if test "x$prefix" = xNONE; then
           ac_bzlib_path="$ac_default_prefix"
        else
           ac_bzlib_path="$prefix"
        fi
    else
        want_bzlib="yes"
        ac_bzlib_path="$withval"
    fi
    ],
    [
    want_bzlib="yes"
    if test "$THIRDPARTY" = "build"; then
        want_bzlib="yes"
        build_bzlib="yes"
        if test "x$prefix" = xNONE; then
           ac_bzlib_path="$ac_default_prefix"
        else
           ac_bzlib_path="$prefix"
        fi
    fi
    ])

if test "x$want_bzlib" = "xyes"; then
        AC_MSG_CHECKING([for BZLIB library])
        succeeded=no
        if test "$ac_bzlib_path" != ""; then
            BZLIB_LDFLAGS="-L$ac_bzlib_path/lib -Wl,-rpath -Wl,$ac_bzlib_path/lib"
            BZLIB_CPPFLAGS="-I$ac_bzlib_path/include"
	    BZLIB_TARGET="$ac_bzlib_path/include/bzlib.h"
            succeeded=yes
        else
            for ac_bzlib_path_tmp in $prefix $THIRDPARTY /groups/algorithm/commonlib/$commonlib_version/$compiler_name /usr /usr/local /opt /opt/local /sw ; do
                  if test -e "$ac_bzlib_path_tmp/include/bzlib.h" && test -r "$ac_bzlib_path_tmp/include/bzlib.h"; then
		      BZLIB_LDFLAGS="-L$ac_bzlib_path_tmp/lib -Wl,-rpath -Wl,$ac_bzlib_path_tmp/lib"
                      BZLIB_CPPFLAGS="-I$ac_bzlib_path_tmp/include"
      	              BZLIB_TARGET="$ac_bzlib_path_tmp/include/bzlib.h"
                      succeeded=yes
                      break;
                  fi
            done
        fi

        if test "$succeeded" != "yes" ; then
                AC_MSG_RESULT([no])
        else
                AC_MSG_RESULT([yes])
                AC_SUBST(BZLIB_CPPFLAGS)
                AC_SUBST(BZLIB_TARGET)
                AC_SUBST(BZLIB_LDFLAGS)
                AC_DEFINE(HAVE_BZLIB,,[Defined if we have BZLIB])
                have_bzlib="yes"
        fi
fi
AM_CONDITIONAL([HAVE_BZLIB], [test "$have_bzlib" = "yes"])
AM_CONDITIONAL([BUILD_BZLIB], [test "$build_bzlib" = "yes"])
])
