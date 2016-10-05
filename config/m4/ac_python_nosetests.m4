#####
#
# SYNOPSIS
#
#   AC_PROG_NOSETESTS()
#
# DESCRIPTION
#
#   This macro searches for a Python Nosetests installation on your system. If
#   found you should call nosetests via $(NOSETESTS). 
#
#   In configure.in, use as:
#
#     AC_PROG_NOSETESTS
#

AC_DEFUN([AC_PROG_NOSETESTS],[
        have_nosetests=no
        AC_ARG_VAR([NOSETESTS], [Override the Nosetests executable to use,
           the default is just 'nosetests'])
        if test "$NOSETESTS" == ""; then
           NOSETESTS=nosetests
        fi
        AC_PATH_PROG([NOSETESTS],[$NOSETESTS])
        if test -z "$NOSETESTS" ; then
           AC_MSG_WARN([cannot find 'nosetests' program.])
           NOSETESTS='echo "Error: nosetests is not installed. " ; false'
        else
          have_nosetests=yes
        fi
AM_CONDITIONAL([HAVE_NOSETESTS], [test "$have_nosetests" = "yes"])
])
