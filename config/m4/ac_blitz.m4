# This looks for the Blitz++ library. If we find them, we set the Makefile
# conditional HAVE_BLITZ. We as set BLITZ_CPPFLAGS and BLITZ_LDFLAGS
# 
# To allow users to build there own copy of Blitz++, we also define
# BUILD_BLITZ

AC_DEFUN([AC_BLITZ],
[
have_blitz="no"
build_blitz="no"
AC_ARG_WITH([blitz],
        AS_HELP_STRING([--with-blitz@<:@=DIR@:>@], [use Blitz++ (default is yes if found) - it is possible to specify the root directory for Blitz++ (optional). You can also specify "build" if you want to build your own local copy. See also THIRDPARTY variable described below.]),
        [
    if test "$withval" = "no"; then
        want_blitz="no"
    elif test "$withval" = "yes"; then
        want_blitz="yes"
        ac_blitz_path=""
    elif test "$withval" = "build"; then
        want_blitz="yes"
        build_blitz="yes"
        if test "x$prefix" = xNONE; then
           ac_blitz_path="$ac_default_prefix"
        else
           ac_blitz_path="$prefix"
        fi
    else
        want_blitz="yes"
        ac_blitz_path="$withval"
    fi
    ],
    [
    want_blitz="yes"
    if test "$THIRDPARTY" = "build"; then
        want_blitz="yes"
        build_blitz="yes"
        if test "x$prefix" = xNONE; then
           ac_blitz_path="$ac_default_prefix"
        else
           ac_blitz_path="$prefix"
        fi
    fi
    ])

if test "x$want_blitz" = "xyes"; then
        AC_MSG_CHECKING([for Blitz library])
        succeeded=no
        if test "$ac_blitz_path" != ""; then
            BLITZ_LDFLAGS="-L$ac_blitz_path/lib -lblitz"
            BLITZ_CPPFLAGS="-I$ac_blitz_path/include"
            succeeded=yes
        else
            for ac_blitz_path_tmp in $prefix $THIRDPARTY /groups/algorithm/commonlib/$commonlib_version/$compiler_name /usr /usr/local /opt /opt/local /sw ; do
                  if test -e "$ac_blitz_path_tmp/include/blitz/array.h" && test -r "$ac_blitz_path_tmp/include/blitz/array.h"; then
                      BLITZ_LDFLAGS="-L$ac_blitz_path_tmp/lib -lblitz"
                      BLITZ_CPPFLAGS="-I$ac_blitz_path_tmp/include"
                      succeeded=yes
                      break;
                  fi
            done
        fi

        if test "$succeeded" != "yes" ; then
                AC_MSG_RESULT([no])
        else
                AC_MSG_RESULT([yes])
                AC_SUBST(BLITZ_CPPFLAGS)
                AC_SUBST(BLITZ_LDFLAGS)
                have_blitz="yes"
        fi
fi
AM_CONDITIONAL([BUILD_BLITZ], [test "$build_blitz" = "yes"])
])