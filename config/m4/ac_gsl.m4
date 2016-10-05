# This looks for the HDF 5 libraries. If we find them, we set the Makefile
# conditional HAVE_GSL. We as set GSL_CPPFLAGS and GSL_LDFLAGS
# 
# To allow users to build there own copy of GSL, we also define
# BUILD_GSL

AC_DEFUN([AC_GSL],
[
have_gsl="no"
build_gsl="no"
AC_ARG_WITH([gsl],
        AS_HELP_STRING([--with-gsl@<:@=DIR@:>@], [use GSL (default is yes if found) - it is possible to specify the root directory for GSL (optional). You can also specify "build" if you want to build your own local copy. See also THIRDPARTY variable described below.]),
        [
    if test "$withval" = "no"; then
        want_gsl="no"
    elif test "$withval" = "yes"; then
        want_gsl="yes"
        ac_gsl_path=""
    elif test "$withval" = "build"; then
        want_gsl="yes"
        build_gsl="yes"
        if test "x$prefix" = xNONE; then
           ac_gsl_path="$ac_default_prefix"
        else
           ac_gsl_path="$prefix"
        fi
    else
        want_gsl="yes"
        ac_gsl_path="$withval"
    fi
    ],
    [
    want_gsl="yes"
    if test "$THIRDPARTY" = "build"; then
        want_gsl="yes"
        build_gsl="yes"
        if test "x$prefix" = xNONE; then
           ac_gsl_path="$ac_default_prefix"
        else
           ac_gsl_path="$prefix"
        fi
    fi
    ])

if test "x$want_gsl" = "xyes"; then
        AC_MSG_CHECKING([for GSL library])
        succeeded=no
        if test "$ac_gsl_path" != ""; then
            GSL_LDFLAGS="-L$ac_gsl_path/lib -lgsl -lgslcblas"
            GSL_CPPFLAGS="-I$ac_gsl_path/include"
            succeeded=yes
        else
            for ac_gsl_path_tmp in $prefix $THIRDPARTY /groups/algorithm/commonlib/$commonlib_version/$compiler_name /usr /usr/local /opt /opt/local /sw ; do
                  if test -e "$ac_gsl_path_tmp/include/gsl/gsl_math.h" && test -r "$ac_gsl_path_tmp/include/gsl/gsl_math.h"; then
                      GSL_LDFLAGS="-L$ac_gsl_path_tmp/lib -lgsl -lgslcblas"
                      GSL_CPPFLAGS="-I$ac_gsl_path_tmp/include"
                      succeeded=yes
                      break;
                  fi
            done
        fi

        if test "$succeeded" != "yes" ; then
                AC_MSG_RESULT([no])
        else
                AC_MSG_RESULT([yes])
                AC_SUBST(GSL_CPPFLAGS)
                AC_SUBST(GSL_LDFLAGS)
                AC_DEFINE(HAVE_GSL,,[Defined if we have GSL])
                have_gsl="yes"
        fi
fi
AM_CONDITIONAL([HAVE_GSL], [test "$have_gsl" = "yes"])
AM_CONDITIONAL([BUILD_GSL], [test "$build_gsl" = "yes"])
])
