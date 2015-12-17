# This looks for the 2stream library. If we find them, we set the Makefile
# conditional HAVE_TWOSTREAM. We also set TWOSTREAM_CPPFLAGS and TWOSTREAM_LDFLAGS
# 

AC_DEFUN([AC_TWOSTREAM],
[
have_twostream="no"
build_twostream="no"
AC_ARG_WITH([2stream],
        AS_HELP_STRING([--with-2stream@<:@=DIR@:>@], [use 2stream (default is yes if found) - it is possible to specify the root directory for 2stream (optional). You can also specify "build" if you want to build your own local copy. See also THIRDPARTY variable described below.]),
        [
    if test "$withval" = "no"; then
        want_twostream="no"
    elif test "$withval" = "yes"; then
        want_twostream="yes"
        ac_twostream_path=""
    elif test "$withval" = "build"; then
        want_twostream="yes"
        build_twostream="yes"
        if test "x$prefix" = xNONE; then
           ac_twostream_path="$ac_default_prefix"
        else
           ac_twostream_path="$prefix"
        fi
    else
        want_twostream="yes"
        ac_twostream_path="$withval"
    fi
    ],
    [
    want_twostream="yes"
    if test "$THIRDPARTY" = "build"; then
        want_twostream="yes"
        build_twostream="yes"
        if test "x$prefix" = xNONE; then
           ac_twostream_path="$ac_default_prefix"
        else
           ac_twostream_path="$prefix"
        fi
    fi
    ])

if test "x$want_twostream" = "xyes"; then
        AC_MSG_CHECKING([for 2stream library])
        succeeded=no
        if test "$ac_twostream_path" != ""; then
            TWOSTREAM_LDFLAGS="$ac_twostream_path/lib/libtwostream.la"
            TWOSTREAM_CPPFLAGS="-I$ac_twostream_path/include"
            succeeded=yes
        else
            for ac_twostream_path_tmp in $prefix $THIRDPARTY /groups/algorithm/commonlib/$commonlib_version/$compiler_name /usr /usr/local /opt /opt/local /sw ; do
                  if test -e "$ac_twostream_path_tmp/lib/libtwostream.la" && test -r "$ac_twostream_path_tmp/lib/libtwostream.la"; then
                      TWOSTREAM_LDFLAGS="$ac_twostream_path_tmp/lib/libtwostream.la"
                      TWOSTREAM_CPPFLAGS="-I$ac_twostream_path_tmp/include"
                      succeeded=yes
                      break;
                  fi
            done
        fi

        if test "$succeeded" != "yes" ; then
                AC_MSG_RESULT([no])
        else
                AC_MSG_RESULT([yes])
                AC_SUBST(TWOSTREAM_CPPFLAGS)
                AC_SUBST(TWOSTREAM_LDFLAGS)
                have_twostream="yes"
        fi
fi
AM_CONDITIONAL([BUILD_TWOSTREAM], [test "$build_twostream" = "yes"])
])
