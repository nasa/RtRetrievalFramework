#####
#
# SYNOPSIS
#
#   AC_PARALLEL()
#
# DESCRIPTION
#
#   This macro searches for a GNU parallel installation on your system. If
#   found you should call parallel via $(PARALLEL). 
#
#   In configure.in, use as:
#
#     AC_PROG_PARALLEL
#

AC_DEFUN([AC_PROG_PARALLEL],[
have_parallel="no"
build_parallel="no"
AC_ARG_WITH([parallel],
        AS_HELP_STRING([--with-parallel@<:@=DIR@:>@], [use PARALLEL (default is yes if found) - it is possible to specify the root directory for PARALLEL (optional). You can also specify "build" if you want to build your own local copy. See also THIRDPARTY variable described below.]),
        [
    if test "$withval" = "no"; then
        want_parallel="no"
    elif test "$withval" = "yes"; then
        want_parallel="yes"
        ac_parallel_path=""
    elif test "$withval" = "build"; then
        want_parallel="yes"
        build_parallel="yes"
        if test "x$prefix" = xNONE; then
           ac_parallel_path="$ac_default_prefix"
        else
           ac_parallel_path="$prefix"
        fi
    else
        want_parallel="yes"
        ac_parallel_path="$withval"
    fi
    ],
    [
    want_parallel="yes"
    if test "$THIRDPARTY" = "build"; then
        want_parallel="yes"
        build_parallel="yes"
        if test "x$prefix" = xNONE; then
           ac_parallel_path="$ac_default_prefix"
        else
           ac_parallel_path="$prefix"
        fi
    fi
    ])
if test "x$want_parallel" = "xyes"; then
   if test "$ac_parallel_path" != ""; then
      PARALLEL=$ac_parallel_path/bin/parallel
      succeeded=yes
   else
      AC_MSG_CHECKING([for parallel])
      succeeded=no
      AC_PATH_PROG([PARALLEL],[parallel])
      if test -z "$PARALLEL" ; then
          have_parallel=no
          PARALLEL='echo "Error: parallel is not installed. " ; false'
      else
         have_parallel=yes
      fi
      if test "$succeeded" != "yes" ; then
           AC_MSG_RESULT([no])
      else
           AC_MSG_RESULT([yes])
      fi
   fi   
   AC_SUBST(PARALLEL)
fi
AM_CONDITIONAL([HAVE_PARALLEL], [test "$have_parallel" = "yes"])
AM_CONDITIONAL([BUILD_PARALLEL], [test "$build_parallel" = "yes"])
])
