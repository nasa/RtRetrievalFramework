#####
#
# SYNOPSIS
#
#   AC_PROG_SPHINX()
#
# DESCRIPTION
#
#   This macro searches for a Python Sphinx installation on your system. If
#   found you should call sphinx-build via $(SPHINXBUILD). 
#
#   In configure.in, use as:
#
#     AC_PROG_SPHINX
#

AC_DEFUN([AC_PROG_SPHINX],[
        have_sphinx=no
        AC_ARG_VAR([SPHINXBUILD], [Override the Sphinx build executable to use,
           the default is just 'sphinx-build'])
        if test "$SPHINXBUILD" == ""; then
           SPHINXBUILD=sphinx-build
        fi
        AC_PATH_PROG([SPHINXBUILD],[$SPHINXBUILD])
        if test -z "$SPHINXBUILD" ; then
           AC_MSG_WARN([cannot find 'sphinx-build' program.])
           SPHINXBUILD='echo "Error: sphinx-build is not installed. " ; false'
        else
          have_sphinx=yes
        fi
AM_CONDITIONAL([HAVE_SPHINX], [test "$have_sphinx" = "yes"])
])
