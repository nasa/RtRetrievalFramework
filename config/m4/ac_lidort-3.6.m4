# This looks for the lidort-3.6 library. If we find them, we set the Makefile
# conditional HAVE_LIDORT_3_6. We also set LIDORT_3_6_CPPFLAGS and LIDORT_3_6_LDFLAGS
# 
# To allow users to build there own copy of LIDORT-3.6, we also define
# BUILD_LIDORT_3_6

AC_DEFUN([AC_LIDORT_3_6],
[
have_lidort_3_6="no"
build_lidort_3_6="no"
AC_ARG_WITH([lidort-3.6],
        AS_HELP_STRING([--with-lidort-3.6@<:@=DIR@:>@], [use LIDORT 3.6 (default is yes if found) - it is possible to specify the root directory for LIDORT 3.6 (optional). You can also specify "build" if you want to build your own local copy. See also THIRDPARTY variable described below.]),
        [
    if test "$withval" = "no"; then
        want_lidort_3_6="no"
    elif test "$withval" = "yes"; then
        want_lidort_3_6="yes"
        ac_lidort_3_6_path=""
    elif test "$withval" = "build"; then
        want_lidort_3_6="yes"
        build_lidort_3_6="yes"
        if test "x$prefix" = xNONE; then
           ac_lidort_3_6_path="$ac_default_prefix"
        else
           ac_lidort_3_6_path="$prefix"
        fi
    else
        want_lidort_3_6="yes"
        ac_lidort_3_6_path="$withval"
    fi
    ],
    [
    want_lidort_3_6="yes"
    if test "$THIRDPARTY" = "build"; then
        want_lidort_3_6="yes"
        build_lidort_3_6="yes"
        if test "x$prefix" = xNONE; then
           ac_lidort_3_6_path="$ac_default_prefix"
        else
           ac_lidort_3_6_path="$prefix"
        fi
    fi
    ])

if test "x$want_lidort_3_6" = "xyes"; then
        AC_MSG_CHECKING([for LIDORT 3.6 library])
        succeeded=no
        if test "$ac_lidort_3_6_path" != ""; then
            LIDORT_3_6_LDFLAGS="$ac_lidort_3_6_path/lib/liblidort-3.6.0.la"
            LIDORT_3_6_CPPFLAGS="-I$ac_lidort_3_6_path/include"
            succeeded=yes
        else
            for ac_lidort_3_6_path_tmp in $prefix $THIRDPARTY /groups/algorithm/commonlib/$commonlib_version/$compiler_name /usr /usr/local /opt /opt/local /sw ; do
                  if test -e "$ac_lidort_3_6_path_tmp/lib/liblidort-3.6.0.la" && test -r "$ac_lidort_3_6_path_tmp/lib/liblidort-3.6.0.la"; then
                      LIDORT_3_6_LDFLAGS="$ac_lidort_3_6_path_tmp/lib/liblidort-3.6.0.la"
                      LIDORT_3_6_CPPFLAGS="-I$ac_lidort_3_6_path_tmp/include"
                      succeeded=yes
                      break;
                  fi
            done
        fi

        if test "$succeeded" != "yes" ; then
                AC_MSG_RESULT([no])
        else
                AC_MSG_RESULT([yes])
                AC_SUBST(LIDORT_3_6_CPPFLAGS)
                AC_SUBST(LIDORT_3_6_LDFLAGS)
                have_lidort_3_6="yes"
        fi
fi
AM_CONDITIONAL([BUILD_LIDORT_3_6], [test "$build_lidort_3_6" = "yes"])
])
