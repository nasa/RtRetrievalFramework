# This looks for the HDF 5 libraries. If we find them, we set the Makefile
# conditional HAVE_HDF5. We as set HDF5_CPPFLAGS and HDF5_LDFLAGS
# 
# To allow users to build there own copy of HDF5, we also define
# BUILD_HDF5

AC_DEFUN([AC_HDF5],
[
have_hdf5="no"
build_hdf5="no"
AC_ARG_WITH([hdf5],
        AS_HELP_STRING([--with-hdf5@<:@=DIR@:>@], [use HDF 5 (default is yes if found) - it is possible to specify the root directory for HDF 5 (optional). You can also specify "build" if you want to build your own local copy. See also THIRDPARTY variable described below.]),
        [
    if test "$withval" = "no"; then
        want_hdf5="no"
    elif test "$withval" = "yes"; then
        want_hdf5="yes"
        ac_hdf5_path=""
    elif test "$withval" = "build"; then
        want_hdf5="yes"
        build_hdf5="yes"
        if test "x$prefix" = xNONE; then
           ac_hdf5_path="$ac_default_prefix"
        else
           ac_hdf5_path="$prefix"
        fi
    else
        want_hdf5="yes"
        ac_hdf5_path="$withval"
    fi
    ],
    [
    want_hdf5="yes"
    if test "$THIRDPARTY" = "build"; then
        want_hdf5="yes"
        build_hdf5="yes"
        if test "x$prefix" = xNONE; then
           ac_hdf5_path="$ac_default_prefix"
        else
           ac_hdf5_path="$prefix"
        fi
    fi
    ])

if test "x$want_hdf5" = "xyes"; then
        AC_MSG_CHECKING([for HDF5 library])
        succeeded=no
        if test "$ac_hdf5_path" != ""; then
            HDF5HOME=$ac_hdf5_path
	    if test -e "$ac_hdf5_path/lib/libhdf5_cpp.la" && test -r "$ac_hdf5_path/lib/libhdf5_cpp.la"; then
	         HDF5_LDFLAGS="$ac_hdf5_path/lib/libhdf5_cpp.la $ac_hdf5_path/lib/libhdf5.la $ac_hdf5_path/lib/libhdf5_hl.la $ac_hdf5_path/lib/libhdf5_hl_cpp.la -L$ac_hdf5_path/lib -lz"
	    else
	        HDF5_LDFLAGS="-L$ac_hdf5_path/lib/lib -lhdf5_cpp -lhdf5 -lhdf5_hl -lhdf5_hl_cpp -lz"
            fi
            HDF5_CPPFLAGS="-I$ac_hdf5_path/include"
            succeeded=yes
        else
            for ac_hdf5_path_tmp in $prefix $THIRDPARTY /groups/algorithm/commonlib/$commonlib_version/$compiler_name /usr /usr/local /opt /opt/local /sw ; do
                  if test -e "$ac_hdf5_path_tmp/lib/libhdf5_cpp.la" && test -r "$ac_hdf5_path_tmp/lib/libhdf5_cpp.la"; then
                      HDF5HOME=$ac_hdf5_path_tmp
                      HDF5_LDFLAGS="-L$ac_hdf5_path_tmp/lib -lhdf5_cpp -lhdf5 -lhdf5_hl -lhdf5_hl_cpp -lz"
                      HDF5_CPPFLAGS="-I$ac_hdf5_path_tmp/include"
                      succeeded=yes
                      break;
                  elif test -e "$ac_hdf5_path_tmp/lib/libhdf5_cpp.so" && test -r "$ac_hdf5_path_tmp/lib/libhdf5_cpp.so"; then
                      HDF5HOME=$ac_hdf5_path_tmp
                      HDF5_LDFLAGS="-L$ac_hdf5_path_tmp/lib -lhdf5_cpp -lhdf5 -lhdf5_hl -lhdf5_hl_cpp -lz"
                      HDF5_CPPFLAGS="-I$ac_hdf5_path_tmp/include"
                      succeeded=yes
                      break;
                  fi
            done
        fi

        if test "$succeeded" != "yes" ; then
                AC_MSG_RESULT([no])
        else
                HDF5_BIN=$HDF5HOME/bin
                AC_MSG_RESULT([yes])
                AC_SUBST(HDF5_CPPFLAGS)
                AC_SUBST(HDF5_LDFLAGS)
                AC_SUBST(HDF5HOME)
                AC_SUBST(HDF5_BIN)
                AC_DEFINE(HAVE_HDF5,,[Defined if we have HDF5])
                have_hdf5="yes"
        fi
fi
if test "x$build_hdf5" = "xyes"; then
   if test "$use_ifort" != "yes" && (test "$use_gfortran" != "yes"); then
      AC_MSG_ERROR([
Right now, HDF5 can only be built with the Intel fortran compiler 
or gfortran. Your compiler isn't one of these. You can try building 
without HDF5.])
   fi
fi
AM_CONDITIONAL([HAVE_HDF5], [test "$have_hdf5" = "yes"])
AM_CONDITIONAL([BUILD_HDF5], [test "$build_hdf5" = "yes"])
])
