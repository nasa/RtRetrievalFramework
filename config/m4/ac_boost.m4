# This looks for the BOOST libraries. If we find them, we set the Makefile
# conditional HAVE_BOOST. We as set BOOST_CPPFLAGS, and BOOST_LIBDIR.
# We also have HAVE_BOOST_STATIC if we have a static version of the 
# libraries.
# 
# To allow users to build there own copy of BOOST, we also define
# BUILD_BOOST

AC_DEFUN([AC_BOOST],
[
have_boost="no"
build_boost="no"
AC_ARG_WITH([boost],
        AS_HELP_STRING([--with-boost@<:@=DIR@:>@], [use BOOST (default is yes if found) - it is possible to specify the root directory for BOOST (optional). You can also specify "build" if you want to build your own local copy. See also THIRDPARTY variable described below.]),
        [
    if test "$withval" = "no"; then
        want_boost="no"
    elif test "$withval" = "yes"; then
        want_boost="yes"
        ac_boost_path=""
    elif test "$withval" = "build"; then
        want_boost="yes"
        build_boost="yes"
        if test "x$prefix" = xNONE; then
           ac_boost_path="$ac_default_prefix"
        else
           ac_boost_path="$prefix"
        fi
    else
        want_boost="yes"
        ac_boost_path="$withval"
    fi
    ],
    [
    want_boost="yes"
    if test "$THIRDPARTY" = "build"; then
        want_boost="yes"
        build_boost="yes"
        if test "x$prefix" = xNONE; then
           ac_boost_path="$ac_default_prefix"
        else
           ac_boost_path="$prefix"
        fi
    fi
    ])

if test "x$want_boost" = "xyes"; then
        AC_MSG_CHECKING([for BOOST library])
        succeeded=no
        if test "$ac_boost_path" != ""; then
            BOOST_LIBDIR="$ac_boost_path/lib"
            BOOST_CPPFLAGS="-I$ac_boost_path/include"
            succeeded=yes
        else
            for ac_boost_path_tmp in $prefix $THIRDPARTY /groups/algorithm/commonlib/$commonlib_version/$compiler_name /usr /usr/local /opt /opt/local /sw ; do
                  if test -e "$ac_boost_path_tmp/include/boost/smart_ptr.hpp" && test -r "$ac_boost_path_tmp/include/boost/smart_ptr.hpp"; then
		      BOOST_LIBDIR="$ac_boost_path_tmp/lib"
                      BOOST_CPPFLAGS="-I$ac_boost_path_tmp/include"
                      succeeded=yes
                      break;
                  fi
            done
        fi

        if test "$succeeded" != "yes" ; then
                AC_MSG_RESULT([no])
        else
                AC_MSG_RESULT([yes])
	        BOOST_CPPFLAGS="$BOOST_CPPFLAGS -DBOOST_TEST_DYN_LINK"
                AC_SUBST(BOOST_CPPFLAGS)
                AC_SUBST(BOOST_LIBDIR)
                AC_DEFINE(HAVE_BOOST,,[Defined if we have BOOST])
                have_boost="yes"
        fi
fi
AM_CONDITIONAL([HAVE_BOOST_STATIC], [test -e "$BOOST_LIBDIR/libboost_regex.a"])
AM_CONDITIONAL([HAVE_BOOST], [test "$have_boost" = "yes"])
AM_CONDITIONAL([BUILD_BOOST], [test "$build_boost" = "yes"])
])
