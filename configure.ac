dnl Process this file with autoconf to produce a configure script.

AC_INIT(full_physics, 11.0.05)

# Base version for commlib.
[commonlib_version='11.0.00a']

AC_COPYRIGHT(
[Copyright 2020, California Institute of Technology. 
ALL RIGHTS RESERVED. U.S. Government Sponsorship acknowledged.])

AC_CONFIG_AUX_DIR([config])
# The obscure looking tar-pax here sets automake to allow file names longer
# than 99 characters to be included in the dist tar. See
# http://noisebleed.blogetery.com/2010/02/27/tar-file-name-is-too-long-max-99/#howtofixit
AM_INIT_AUTOMAKE([1.9 tar-pax])
AM_MAINTAINER_MODE
AC_PROG_LIBTOOL

#======================================================================
# Build our own version of the python modules we depend on. Useful for
# installing on Pleiades
#======================================================================

AC_ARG_WITH([python-module],
           [AS_HELP_STRING([--with-python-module],
             [Build our own copy of the python modules we depend on. Normally you would just install this modules on your system, but for Pleiades this is convenient to build ourselves..])],
           [build_nose=no
	    build_h5py=yes
	    build_future=no
	    build_six=no
	    build_python=yes
	    build_ply=yes
            build_ruport=yes
	    rubygemdir='${prefix}'`${RUBY} -rrbconfig -rrubygems -e 'print Gem.dir.sub(RbConfig::CONFIG.fetch(%q(prefix)), %q())'`
	   ],
           [build_nose=no
	    build_h5py=no
	    build_six=no
	    build_future=no
	    build_python=no
	    build_ply=yes
	    build_ply=no
            build_ruport=no
            rubygemdir=""
	   ])
AM_CONDITIONAL([BUILD_NOSE], [test x$build_nose = xyes])
AM_CONDITIONAL([BUILD_H5PY], [test x$build_h5py = xyes])
AM_CONDITIONAL([BUILD_PYTHON], [test x$build_python = xyes])
AM_CONDITIONAL([BUILD_FUTURE], [test x$build_future = xyes])
AM_CONDITIONAL([BUILD_SIX], [test x$build_six = xyes])
AM_CONDITIONAL([BUILD_PLY], [test x$build_ply = xyes])
AM_CONDITIONAL([BUILD_RUPORT], [test x$build_ruport = xyes])
AC_SUBST(rubygemdir)

AC_PREFIX_DEFAULT([`pwd`/install])
AC_PROG_CC
AC_PROG_FC([ifort gfortran f90 f95 g95 pathf90 pgf90])
AC_PROG_F77([${FC} ifort gfortran f90 f95 g95 pathf90 pgf90])
AC_PROG_FC_C_O
AC_PROG_CXX
AC_FC_LIBRARY_LDFLAGS
AC_DETERMINE_COMPILER
AC_PROG_SPHINX
AC_PROG_PARALLEL
AC_PROG_NOSETESTS
if test "$have_nosetests" = "no"; then
   AC_MSG_WARN([
The python nosetests package is not found. This is used to run some unit tests.
Without this package, these unit tests won't be run. Try specifying
location using NOSETESTS variable if configure couldn't find the
package. 
])
fi
# Latest swig as of 11/15 is 3.0.7. Go ahead and require this, just so we
# aren't dealing with any backward bugs etc. Since this is easy enough to
# install, no reason not to require the latest
AC_PROG_SWIG(3.0.7)
SWIG_ENABLE_CXX
# 2.6 is required for SWIG. Note that we check again for any version of
# python below, we can build the rest of the system with an older version
# of python
AC_PYTHON_DEVEL(>= '2.6')
SWIG_PYTHON
AC_RUBY
if test $PYTHON == ":"; then
   AC_MSG_ERROR([
We currently require python to build the Full Physics library.
Try specifying the location using --with-python if configure
couldn't find python. If you don't have it on your system, try
installing it.])
fi

if test $RUBY == ":"; then
   AC_MSG_ERROR([
We currently require ruby to build the Full Physics library.
Try specifying the location using RUBY=<location> if configure
couldn't find ruby. If you don't have it on your system, try
installing it.])
fi

AC_COMPILER_FLAG
AC_ARG_VAR([THIRDPARTY], [
It is common that either 1) you want to build a local copy of all
thirdparty software or 2) you have a central location not already
searched for that you want to be searched. You can specify

THIRDPARTY=build or THIRDPARTY=<directory>. This gets added to all
the third party library checks.
])
AC_HDF5
AC_LIDORT_3_8
if test "$have_lidort_3_8" = "no"; then
   AC_MSG_ERROR([
The LIDORT 3.8 library is required by Full Physics. Try specifying
location using --with-lidort-3.8 if configure couldn't find the
library. You can specify --with-lidort-3.8=build if you want to build your
own local copy of LIDORT 3.8.

You can also specify THIRDPARTY=build if you want to build a local
copy of all the Thirdparty libraries.
])
fi
AC_TWOSTREAM
if test "$have_twostream" = "no"; then
   AC_MSG_ERROR([
The 2stream library is required by Full Physics. Try specifying
location using --with-2stream if configure couldn't find the
library. You can specify --with-2stream=build if you want to build your
own local copy of 2stream.

You can also specify THIRDPARTY=build if you want to build a local
copy of all the Thirdparty libraries.
])
fi
AC_BLITZ
if test "$have_blitz" = "no"; then
   AC_MSG_ERROR([
The Blitz library is required by Full Physics. Try specifying
location using --with-blitz if configure couldn't find the
library. You can specify --with-blitz=build if you want to build your
own local copy of Blitz++.

You can also specify THIRDPARTY=build if you want to build a local
copy of all the Thirdparty libraries.
])
fi
AC_BOOST
if test "$have_boost" = "no"; then
   AC_MSG_ERROR([
The BOOST library is required by Full Physics. Try specifying
location using --with-boost if configure couldn't find the
library. You can specify --with-boost=build if you want to build your
own local copy of BOOST.

You can also specify THIRDPARTY=build if you want to build a local
copy of all the Thirdparty libraries.
])
fi
AC_BZLIB
if test "$build_boost" = "yes"; then
   if test "$have_bzlib" = "no"; then
      AC_MSG_ERROR([
The bzlib library is required by Full Physics if you are building the boost
library. Try specifying location using --with-bzlib if configure couldn't 
find the library. You can specify --with-bzlib=build if you want to build your
own local copy of bzlib.

You can also specify THIRDPARTY=build if you want to build a local
copy of all the Thirdparty libraries.
])
   fi
fi
AC_GSL
if test "$have_gsl" = "no"; then
   AC_MSG_ERROR([
The GSL library is required by Full Physics. Try specifying
location using --with-gsl if configure couldn't find the
library. You can specify --with-gsl=build if you want to build your
own local copy of GSL.

You can also specify THIRDPARTY=build if you want to build a local
copy of all the Thirdparty libraries.
])
fi

AC_LUA

# We don't actually use any threading code, but on some systems the static
# libraries such as libgfortran require pthread library. So go ahead and
# check for these libraries, and if found we include in our linking.
ACX_PTHREAD

DX_INIT_DOXYGEN(FullPhysics, config/doxygen.cfg, doc)

# check if we have feenableexcept available. This is required by C99, but
# isn't on the MAC for some reason.
save_LIBS="$LIBS"; LIBS="$LIBS -lm"
AC_CHECK_FUNCS([feenableexcept], AC_DEFINE(HAVE_FEENABLEEXCEPT,1,[Define if you have feenableexcept]))
LIBS="$save_LIBS"

#======================================================================
# We use Fortran 2003, and we know certain compilers don't work. Give
# a warning message if the user tries to use them.
#======================================================================

if test "$use_absoft" = "yes"; then
  AC_MSG_ERROR([
The Full Physics code makes use of a number of Fortran 2003 features.
The Absoft compiler doesn't support these feature. You'll need to use
a different compiler (e.g., ifort or gfortran)
])
fi

#======================================================================
# There a some minor differences in building on a Mac vs. Linux. Most
# of the Makefile doesn't care about this, but there are a few places we
# need to know if we are on a Mac.
#======================================================================

AM_CONDITIONAL([IS_MAC], [test ! -z "`uname | grep Darwin`"])

#======================================================================
# There are a number of directories that we need to keep track of when
# building. It is somewhat arbitrary how we do this. We have choosen
# to collect all of this in the top level configure.in file. We could
# have just as easily specified this is each of our automake files,
# but we needed to pick a way.
#======================================================================

if test "x$prefix" = xNONE; then
   prefix_value="$ac_default_prefix"
else
   prefix_value="$prefix"
fi

AC_SUBST([docpythonsrc], [doc/python])
AC_SUBST([pythonsrc], [bindings/python])
AC_SUBST([swigsrc], [bindings/python/swig])
AC_SUBST([srcscript], [script])
AC_SUBST([binsrc], [bin])
AC_SUBST([thirdpartysrc], [thirdparty])
AC_SUBST([lidort_3_8_src], [thirdparty/lidort-3.8])
AC_SUBST([twostream_src], [thirdparty/2stream])
AC_SUBST([supportsrc], [lib/Support])
AC_SUBST([interfacesrc], [lib/Interface])
AC_SUBST([implsrc], [lib/Implementation])
AC_SUBST([fixedlevelsrc], [lib/FixedLevel])
AC_SUBST([deprecatedsrc], [lib/Deprecated])
AC_SUBST([luasrc], [lib/Lua])
AC_SUBST([outsrc], [lib/RegisterOutput])
AC_SUBST([pythonlibsrc], [lib/Python])
AC_SUBST([pythontemplatesrc], [lib/Python/template])
AC_SUBST([l2analysissrc], [lib/Python/l2_analysis])
AC_SUBST([populateconfigsrc], [lib/Python/populate_config])
AC_SUBST([supportutilssrc], [support/utils])
AC_SUBST([fpsrc], [lib])
AC_SUBST([fileiosrc], [lib/fileio])
AC_SUBST([l_radsrc], [lib/Implementation/l_rad_fortran])
AC_SUBST([full_run], [unit_test_data/full_run])
AC_SUBST([fts_run], [unit_test_data/fts_run])
AC_SUBST([configdir], [${prefix_value}/etc/full_physics/config])
AC_SUBST([inputdir], [${prefix_value}/etc/full_physics/input])
AC_SUBST([pkgpythontemplatedir], [${pkgpythondir}/template])
AC_SUBST([commonconfigsrc], [input/common/config])
AC_SUBST([commoninputsrc], [input/common/input])
AC_SUBST([ftsconfigsrc], [input/fts/config])
AC_SUBST([ftsinputsrc], [input/fts/input])
AC_SUBST([gosatconfigsrc], [input/gosat/config])
AC_SUBST([gosatinputsrc], [input/gosat/input])
AC_SUBST([ococonfigsrc], [input/oco/config])
AC_SUBST([ocoinputsrc], [input/oco/input])


#======================================================================
# Option to bypass installing documentation. This is useful
# for a build environment that is used for development, but you'll want all
# of your "real" build environments to also install the documentation.
#======================================================================

AC_ARG_WITH([documentation],
           [AS_HELP_STRING([--without-documentation],
             [Bypass installing documentation. This is useful for a build environment that is used for development, but you'll want all of your "real" build 
environments to also install the documentation.])],
           [],
           [with_documentation=no])
AM_CONDITIONAL([WITH_DOCUMENTATION], [test x$with_documentation = xyes])

#======================================================================
# Allow location of absco tables to be changed. This is used by the
# l2_fp_run test, as well as some unit tests.
#======================================================================

AC_ARG_WITH([absco],
AS_HELP_STRING([--with-absco@<:@=DIR@:>@], [give the location of the ABSCO tables to use in the l2_fp_run and unit tests. The default is /groups/algorithm/l2_fp/absco. This is only needed for the l2_fp_run and some unit tests, if you don't have the ABSCO tables available you can still build everything else, just not do the l2_fp_run and some unit tests.]), [ ac_absco_path="$withval" ], [ ac_absco_path="/groups/algorithm/l2_fp/absco"] )
AC_SUBST([abscodir], ["$ac_absco_path"])
if test ! -d "$ac_absco_path"; then
   AC_MSG_WARN([The ABSCO directory $ac_absco_path is not found. This is only needed for the l2_fp_run and some unit tests, if you don't have the ABSCO tables available you can still build everything else, just not do the l2_fp_run and some unit tests. You can specify an ABSCO directory if desired using --with-absco.])
fi

#======================================================================
# Allow location of merra tables to be changed. This is used by the
# l2_fp_run test, as well as some unit tests.
#======================================================================

AC_ARG_WITH([merra],
AS_HELP_STRING([--with-merra@<:@=DIR@:>@], [give the location of the MERRA tables to use in the l2_fp_run and unit tests. The default is /groups/algorithm/l2_fp/merra_composite. This is only needed for the l2_fp_run and some unit tests, if you don't have the MERRA tables available you can still build everything else, just not do the l2_fp_run and some unit tests.]), [ ac_merra_path="$withval" ], [ ac_merra_path="/groups/algorithm/l2_fp/merra_composite"] )
AC_SUBST([merradir], ["$ac_merra_path"])
if test ! -d "$ac_merra_path"; then
   AC_MSG_WARN([The MERRA directory $ac_merra_path is not found. This is only needed for the l2_fp_run and some unit tests, if you don't have the MERRA tables available you can still build everything else, just not do the l2_fp_run and some unit tests. You can specify an MERRA directory if desired using --with-merra.])
fi

AC_CONFIG_FILES([Makefile script/setup_fp_env.sh lib/Support/global_fixture_default.cc])
AC_CONFIG_FILES([script/python_wrap.sh], [chmod +x script/python_wrap.sh])


if test "$have_swig" = "yes" -a "$USE_MAINTAINER_MODE" != "no"; then
  generate_swig="yes"
else
  generate_swig="no" 
fi

#======================================================================
# Require doxygen if we are generating SWIG.
#======================================================================

if test "$generate_swig" = "yes" -a -z "$DX_DOXYGEN"; then
   AC_MSG_ERROR([doxygen is required if we are building SWIG])
fi


#======================================================================
# Allow max layers in LIDORT the be changed.
#======================================================================

AC_ARG_WITH([lidort-maxlayer],
AS_HELP_STRING([--with-lidort-maxlayer=value], [optionally change the maximum number of layers we can use in LIDORT. The default is 23. Note that this can increase the memory use a great deal if you set this to a large number (which is why we don't always have this set to a large number). The actual number of layers used in a particular run can be smaller than this, this is just the maximum allowed value.]), [ ac_lidort_maxlayer="$withval" ], [ ac_lidort_maxlayer="23"] )
AC_SUBST([lidort_maxlayer], ["$ac_lidort_maxlayer"])

#======================================================================
# Allow max atmoswfs in LIDORT the be changed.
#======================================================================

AC_ARG_WITH([lidort-maxatmoswfs],
AS_HELP_STRING([--with-lidort-maxatmoswfs=value], [optionally change the maximum number of atmospheric weighting functions (jacobians) we can use in LIDORT. The default is 7. Note that this can increase the memory use a great deal if you set this to a large number (which is why we don't always have this set to a large number). The actual number of weigthing functions used in a particular run can be smaller than this, this is just the maximum allowed value.]), [ ac_lidort_maxatmoswfs="$withval" ], [ ac_lidort_maxatmoswfs="7"] )
AC_SUBST([lidort_maxatmoswfs], ["$ac_lidort_maxatmoswfs"])

#======================================================================
# Generate lidort parameters with this maximum value.
# We have extra logic here to update the file only if it doesn't exist
# or is changed, without this we end up regenerating lidort and everything
# that depends on it every time we rerun the configuration script.
#======================================================================

AC_CONFIG_FILES([thirdparty/lidort-3.8/lidort_def/lidort_pars.F90.tmp:thirdparty/lidort-3.8/lidort_def/lidort_pars.F90.in],
[
 f1=thirdparty/lidort-3.8/lidort_def/lidort_pars.F90.tmp
 f2=thirdparty/lidort-3.8/lidort_def/lidort_pars.F90
 AS_IF([test -f $f2],
 [
   AS_IF([cmp -s $f1 $f2 2>/dev/null],
   [
      AC_MSG_NOTICE([$f2 is unchanged])
   ],[
      AC_MSG_NOTICE([creating $f2])
      rm -f $f2
      cp $f1 $f2
   ])
  ],[
    AC_MSG_NOTICE([creating $f2])
    rm -f $f2
    cp $f1 $f2
  ])
])

#======================================================================
# Move CXXFLAGS flags into AM_CXXFLAGS so we can override this if 
# desired
#======================================================================

AM_CXXFLAGS="$CXXFLAGS"
AC_SUBST([AM_CXXFLAGS])
CXXFLAGS=""

AC_OUTPUT

#======================================================================
# Give a nice log message describing what has been configured, and which
# optional libraries are going to be used.
#======================================================================

LOC_MSG()
LOC_MSG([Level 2 Full Physics is now configured])
LOC_MSG()
LOC_MSG([  Installation directory:       ${prefix}])
LOC_MSG([  Build debug version:          ${enable_debug}])
LOC_MSG([  Fortran compiler type:        ${compiler_name}])
LOC_MSG([  Fortran compiler:             ${FC} ${FCFLAGS}])
LOC_MSG([  C compiler:                   ${CC} ${CFLAGS}])
LOC_MSG([  CXX compiler:                 ${CXX} ${AM_CXXFLAGS}])
LOC_MSG([  LD Flags:                     ${LDFLAGS}])
LOC_MSG()
LOC_MSG([  HDF5 support:                 ${have_hdf5}])
LOC_MSG([  Lua support:                  ${have_lua}])
LOC_MSG([  Build own HDF5 library:       ${build_hdf5}])
LOC_MSG([  Build own LIDORT 3.8 library: ${build_lidort_3_8}])
LOC_MSG([  Build own 2stream library:    ${build_twostream}])
LOC_MSG([  Build own BZIB library:       ${build_bzlib}])
LOC_MSG([  Build own BOOST library:      ${build_boost}])
LOC_MSG([  Build own Blitz++ library:    ${build_blitz}])
LOC_MSG([  Build own GSL library:        ${build_gsl}])
LOC_MSG([  Build own LUA library:        ${build_lua}])
LOC_MSG([  Build own GNU parallel:       ${build_parallel}])
LOC_MSG([  Build own nose module:        ${build_nose}])
LOC_MSG([  Build own h5py module:        ${build_h5py}])
LOC_MSG([  Build own future module:      ${build_future}])
LOC_MSG([  Build own six module:         ${build_six}])
LOC_MSG([  Build own PLY module:         ${build_ply}])
LOC_MSG([  Build own Ruport gem:         ${build_ruport}])
LOC_MSG([  Generate SWIG wrappers:       ${generate_swig}])
LOC_MSG([  Compile Python SWIG wrappers: ${build_python_swig}])
LOC_MSG([  Have python Sphinx:           ${have_sphinx}])
LOC_MSG([  Have python nosetests:        ${have_nosetests}])
LOC_MSG([  Install documentation:        ${with_documentation}])

