##### http://autoconf-archive.cryp.to/ac_python_devel.html
#
# SYNOPSIS
#
#   AC_PYTHON_DEVEL([version])
#
# DESCRIPTION
#
#   Note: Defines as a precious variable "PYTHON_VERSION". Don't
#   override it in your configure.ac.
#
#   This macro checks for Python and tries to get the include path to
#   'Python.h'. It provides the $(PYTHON_CPPFLAGS) and
#   $(PYTHON_LDFLAGS) output variables. It also exports
#   $(PYTHON_EXTRA_LIBS) and $(PYTHON_EXTRA_LDFLAGS) for embedding
#   Python in your code.
#
#   You can search for some particular version of Python by passing a
#   parameter to this macro, for example ">= '2.3.1'", or "== '2.4'".
#   Please note that you *have* to pass also an operator along with the
#   version to match, and pay special attention to the single quotes
#   surrounding the version number. Don't use "PYTHON_VERSION" for
#   this: that environment variable is declared as precious and thus
#   reserved for the end-user.
#
#   This macro should work for all versions of Python >= 2.1.0. As an
#   end user, you can disable the check for the python version by
#   setting the PYTHON_NOVERSIONCHECK environment variable to something
#   else than the empty string.
#
#   If you need to use this macro for an older Python version, please
#   contact the authors. We're always open for feedback.
#
# LAST MODIFICATION
#
#   2007-07-31
#
# COPYLEFT
#
#   Copyright (c) 2007 Sebastian Huber <sebastian-huber@web.de>
#   Copyright (c) 2007 Alan W. Irwin <irwin@beluga.phys.uvic.ca>
#   Copyright (c) 2007 Rafael Laboissiere <rafael@laboissiere.net>
#   Copyright (c) 2007 Andrew Collier <colliera@ukzn.ac.za>
#   Copyright (c) 2007 Matteo Settenvini <matteo@member.fsf.org>
#   Copyright (c) 2007 Horst Knorr <hk_classes@knoda.org>
#
#   This program is free software: you can redistribute it and/or
#   modify it under the terms of the GNU General Public License as
#   published by the Free Software Foundation, either version 3 of the
#   License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#   General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see
#   <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright
#   owner gives unlimited permission to copy, distribute and modify the
#   configure scripts that are the output of Autoconf when processing
#   the Macro. You need not follow the terms of the GNU General Public
#   License when using or distributing such scripts, even though
#   portions of the text of the Macro appear in them. The GNU General
#   Public License (GPL) does govern all other use of the material that
#   constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the
#   Autoconf Macro released by the Autoconf Macro Archive. When you
#   make and distribute a modified version of the Autoconf Macro, you
#   may extend this special exception to the GPL to apply to your
#   modified version as well.
#
# This defines the variable build_python_swig and Makefile configure 
# BUILD_PYTHON_SWIG if we find python on the system.

AC_DEFUN([AC_PYTHON_DEVEL],[
build_python_swig="no"

	#
	# Allow the use of a (user set) custom python version
	#
	AC_ARG_VAR([PYTHON],[Override the Python executable to use,
                   the default is just 'python'])
	AC_ARG_VAR([PYTHON_VERSION],[The installed Python
		version to use, for example '2.3'. This string
		will be appended to the Python interpreter
		canonical name.])

AC_ARG_WITH([python-swig],
	AS_HELP_STRING([--with-python-swig], [build SWIG wrappers for use with Python (default is no) - it is possible to specify a path to python (optional)]),
    [
    if test "$withval" = "no"; then
	want_python="no"
    elif test "$withval" = "yes"; then
        want_python="yes"
        ac_python_path=""
    else
        want_python="yes"
        ac_python_path="$withval"
    fi
    ],
    [want_python="no"])

  if test "$build_python" == "yes"; then
     AM_PATH_PYTHON(,, [:])
     # Will need to update this if we change the version we are building
     python_inc_path=python3.9
     python_lib_path=python3.9
     PYTHON=`pwd`"/script/python_wrap.sh" 
     PYTHON_CPPFLAGS="-I\${prefix}/include/${python_inc_path}"
     numpy_cppflags="-I\${prefix}/lib/${python_lib_path}/site-packages/numpy/core/include"
     pythondir="\${prefix}/lib/${python_lib_path}/site-packages" 
     platpythondir="\${prefix}/lib/${python_lib_path}/site-packages" 
     build_python_swig="yes"
     PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $numpy_cppflags"
     AC_SUBST(PYTHON_CPPFLAGS)
     AC_SUBST([platpythondir])
     PYTHON_LDFLAGS=""
     AC_SUBST(PYTHON_LDFLAGS)
   else
     if test "x$want_python" = "xyes"; then
	AC_PATH_PROG([PYTHON],[python[$PYTHON_VERSION]])
	if test "$ac_python_path" != ""; then
	   PYTHON=$ac_python_path
        fi
	if test -z "$PYTHON"; then
	   AC_MSG_RESULT([Cannot find python$PYTHON_VERSION in your system path])
	   PYTHON_VERSION=""
	else
          build_python_swig="yes"

        # If default python isn't new enough, then try again for python2.6
	  if test -n "$1"; then
		AC_MSG_CHECKING([for a version of Python $1])
		ac_supports_python_ver=`$PYTHON -c "import sys; \
			ver = sys.version.split()[[0]]; \
			print (ver $1)"`
		if test "$ac_supports_python_ver" = "False"; then
	           AC_PATH_PROG([PYTHON26],[python2.6])
                   if test ! -z "$PYTHON26"; then
                      PYTHON="$PYTHON26"
                   fi
                fi
          fi


	#
	# Check for a version of Python >= 2.1.0
	#
  	  AC_MSG_CHECKING([for a version of Python >= '2.1.0'])
  	  ac_supports_python_ver=`$PYTHON -c "import sys; \
		ver = sys.version.split()[[0]]; \
		print (ver >= '2.1.0')"`
	  if test "$ac_supports_python_ver" != "True"; then
		if test -z "$PYTHON_NOVERSIONCHECK"; then
			AC_MSG_RESULT([no])
			AC_MSG_WARN([
This version of the AC@&t@_PYTHON_DEVEL macro
doesn't work properly with versions of Python before
2.1.0. You may need to re-run configure, setting the
variables PYTHON_CPPFLAGS, PYTHON_LDFLAGS, PYTHON_SITE_PKG,
PYTHON_EXTRA_LIBS and PYTHON_EXTRA_LDFLAGS by hand.
Moreover, to disable this check, set PYTHON_NOVERSIONCHECK
to something else than an empty string.
])
                        build_python_swig="no"
		else
			AC_MSG_RESULT([skip at user request])
		fi
	  else
		AC_MSG_RESULT([yes])
	  fi

	#
	# if the macro parameter ``version'' is set, honour it
	#
	  if test -n "$1"; then
		AC_MSG_CHECKING([for a version of Python $1])
		ac_supports_python_ver=`$PYTHON -c "import sys; \
			ver = sys.version.split()[[0]]; \
			print (ver $1)"`
		if test "$ac_supports_python_ver" = "True"; then
	   	   AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
			AC_MSG_WARN([this package requires Python $1
to build the python wrappers. The rest of the system can be built
with an older version, so you can ignore this warning if you aren't
using the wrappers. If you have it installed, but it isn't the 
default Python interpreter in your system path, please pass the 
PYTHON_VERSION variable to configure. See ``configure --help'' for reference.
])
			PYTHON_VERSION=""
		fi
	  fi

	#
	# Check if you have distutils, else fail
	#
  	  AC_MSG_CHECKING([for the distutils Python package])
 	  ac_distutils_result=`$PYTHON -c "import distutils" 2>&1`
	  if test -z "$ac_distutils_result"; then
		AC_MSG_RESULT([yes])
	  else
		AC_MSG_RESULT([no])
		AC_MSG_WARN([cannot import Python module "distutils".
Please check your Python installation. The error was:
$ac_distutils_result])
                build_python_swig="no"
		PYTHON_VERSION=""
	  fi

	#
	# Check for Python include path
	#
	  AC_MSG_CHECKING([for Python include path])
	  if test -z "$PYTHON_CPPFLAGS"; then
		python_path=`$PYTHON -c "import distutils.sysconfig; \
           		print (distutils.sysconfig.get_python_inc());"`
		if test -n "${python_path}"; then
		   	python_path="-I$python_path"
		fi
		PYTHON_CPPFLAGS=$python_path
	  fi
	  AC_MSG_RESULT([$PYTHON_CPPFLAGS])
	  AC_SUBST([PYTHON_CPPFLAGS])

          AC_MSG_CHECKING([for numpy include path])
 	  ac_numpy_result=`$PYTHON -c "import numpy.distutils" 2>&1`
	  if test -z "$ac_numpy_result"; then
		AC_MSG_RESULT([yes])
                numpy_cppflags=`$PYTHON -c "import numpy.distutils;  print ('-I' + ' -I'.join(numpy.distutils.misc_util.get_numpy_include_dirs()));"`
                PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $numpy_cppflags"
	  else
		AC_MSG_RESULT([no])
		AC_MSG_WARN([cannot import Python module "numpy.distutils".
Won't be able to build swig wrappers. The error was:
$ac_numpy_result])
                build_python_swig="no"
	  fi
          

	#
	# Check for Python library path
	#
	  AC_MSG_CHECKING([for Python library path])
	  if test -z "$PYTHON_LDFLAGS"; then
		# (makes two attempts to ensure we've got a version number
		# from the interpreter)
		py_version=`$PYTHON -c "from distutils.sysconfig import *; \
			print(' '.join(get_config_vars('VERSION')))"`
		if test "$py_version" == "[None]"; then
			if test -n "$PYTHON_VERSION"; then
				py_version=$PYTHON_VERSION
			else
				py_version=`$PYTHON -c "import sys; \
					print (sys.version[[:3]])"`
			fi
		fi

		PYTHON_LDFLAGS=`$PYTHON -c "from distutils.sysconfig import *; \
			print ('-L/usr/lib -L/usr/lib64 -L' + get_python_lib(0,1) + ' -L' + get_python_lib(0,1) + '/config', \
		      	'-lpython');"`$py_version
	fi
	AC_MSG_RESULT([$PYTHON_LDFLAGS])
	AC_SUBST([PYTHON_LDFLAGS])

	#
	# Check for site packages
	#
	AC_MSG_CHECKING([for Python site-packages path])
	if test -z "$PYTHON_SITE_PKG"; then
		PYTHON_SITE_PKG=`$PYTHON -c "import distutils.sysconfig; \
		        print (distutils.sysconfig.get_python_lib(0,0));"`
	fi
	AC_MSG_RESULT([$PYTHON_SITE_PKG])
	AC_SUBST([PYTHON_SITE_PKG])

	#
	# libraries which must be linked in when embedding
	#
	AC_MSG_CHECKING(python extra libraries)
	if test -z "$PYTHON_EXTRA_LIBS"; then
	   PYTHON_EXTRA_LIBS=`$PYTHON -c "import distutils.sysconfig; \
                conf = distutils.sysconfig.get_config_var; \
                print (conf('LOCALMODLIBS'), conf('LIBS'))"`
	fi
	AC_MSG_RESULT([$PYTHON_EXTRA_LIBS])
	AC_SUBST(PYTHON_EXTRA_LIBS)

	#
	# linking flags needed when embedding
	#
	AC_MSG_CHECKING(python extra linking flags)
	if test -z "$PYTHON_EXTRA_LDFLAGS"; then
		PYTHON_EXTRA_LDFLAGS=`$PYTHON -c "import distutils.sysconfig; \
			conf = distutils.sysconfig.get_config_var; \
			print (conf('LINKFORSHARED'))"`
	fi
	AC_MSG_RESULT([$PYTHON_EXTRA_LDFLAGS])
	AC_SUBST(PYTHON_EXTRA_LDFLAGS)

	#
	# final check to see if everything compiles alright
	#
	AC_MSG_CHECKING([consistency of all components of python development environment])
	AC_LANG_PUSH([C])
	# save current global flags
	LIBS="$ac_save_LIBS $PYTHON_LDFLAGS"
	CPPFLAGS="$ac_save_CPPFLAGS $PYTHON_CPPFLAGS"
        [pythonexists=yes]
# 	AC_TRY_LINK([
# 		#include <Python.h>
# 	],[
# 		Py_Initialize();
# 	],[pythonexists=yes],[pythonexists=no])

	AC_MSG_RESULT([$pythonexists])

        if test ! "$pythonexists" = "yes"; then
	   AC_MSG_WARN([
  Could not link test program to Python. Maybe the main Python library has been
  installed in some non-standard library path. If so, pass it to configure,
  via the LDFLAGS environment variable.
  Example: ./configure LDFLAGS="-L/usr/non-standard-path/python/lib"
  ============================================================================
   You probably have to install the development version of the Python package
   for your distribution.  The exact name of this package varies among them.
  ============================================================================
	   ])
          build_python_swig="no"
	  PYTHON_VERSION=""
	fi
	AC_LANG_POP
	# turn back to default flags
	CPPFLAGS="$ac_save_CPPFLAGS"
	LIBS="$ac_save_LIBS"
	fi
	pythondir=`$PYTHON -c "from distutils.sysconfig import *; print(get_python_lib(False,False,''))"`
	platpythondir=`$PYTHON -c "from distutils.sysconfig import *; print(get_python_lib(True,False,''))"`
	pythondir="\${prefix}/${pythondir}"
	platpythondir="\${prefix}/${platpythondir}"
	AC_SUBST([platpythondir])
      else
        AM_PATH_PYTHON(,, [:])
    fi
fi
AC_SUBST([pkgpythondir], [\${pythondir}/$PACKAGE])
AM_CONDITIONAL([BUILD_PYTHON_SWIG], [test "$build_python_swig" = "yes"])
])
