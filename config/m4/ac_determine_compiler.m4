# The Fortran compiler we are using determines what set of flags we use,
# and in some cases where to look for files (e.g, the HDF 5 libraries are
# separated by compiler type).
#
# This routine files in the automake variables USE_IFORT, USE_ABSOFT, 
# USE_NAG, USE_G95, USE_GFORTRAN, USE_PATHSCALE, USE_PORTLAND, USE_UNKNOWN. 
# One of these flags will be true - the default USE_UNKNOWN is set if we
# don't recongnize the compiler.
#
# We also define the lower case config variables with the same names, e.g.
# $use_ifort, $use_absoft, etc. These are "no" or "yes".
#
# Finally $compiler_name is set to the name of the compiler type.

AC_DEFUN([AC_DETERMINE_COMPILER],
[
use_ifort="no"
use_absoft="no"
use_nag="no"
use_g95="no"
use_gfortran="no"
use_pathscale="no"
use_portland="no"
use_unknown="no"

AC_MSG_CHECKING([compiler type])
# Freeze FC to compiler that is found in the path. This is for compilers
# like ifort which typically have multiple ones on the system. We select
# the compiler at the configure time, and use this one throughout.
FC=`which $FC`
if ($FC -v 2>&1 | grep -q gcc) && (echo $FC | grep -q g95); then
   use_g95="yes"
   compiler_name="g95"
elif ($FC -v 2>&1 | grep -q gcc) && (echo $FC | grep -q gfortran); then
   use_gfortran="yes"
   compiler_name="gfortran"
elif ($FC -V 2>&1 | grep -q "Intel(R) Fortran"); then
   use_ifort="yes"
   compiler_name="ifort"
elif ($FC -V 2>&1 | grep -q "NAGWare Fortran 95"); then
   use_nag="yes"
   compiler_name="NAG"
elif ($FC -V 2>&1 | grep -q "Absoft"); then
   use_absoft="yes"
   compiler_name="Absoft"
elif (echo $FC | grep -q "pathf90"); then
   use_pathscale="yes"
   compiler_name="pathscale"
elif (echo $FC | grep -q "pgf90"); then
   use_pathscale="yes"
   compiler_name="portland group"
else
   use_unknown="yes"
   compiler_name="unknown"
fi
AC_MSG_RESULT([$compiler_name])
AM_CONDITIONAL([USE_IFORT], [test "$use_ifort" = "yes"])
AM_CONDITIONAL([USE_ABSOFT], [test "$use_absoft" = "yes"])
AM_CONDITIONAL([USE_NAG], [test "$use_nag" = "yes"])
AM_CONDITIONAL([USE_G95], [test "$use_g95" = "yes"])
AM_CONDITIONAL([USE_GFORTRAN], [test "$use_gfortran" = "yes"])
AM_CONDITIONAL([USE_PATHSCALE], [test "$use_pathscale" = "yes"])
AM_CONDITIONAL([USE_PORTLAND], [test "$use_portland" = "yes"])
AM_CONDITIONAL([USE_UNKNOWN], [test "$use_unknown" = "yes"])
])


