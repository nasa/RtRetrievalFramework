# This sets up the Fortran compiler flags. You should make sure to
# call AC_DETERMINE_COMPILER before this function. This sets FCFLAGS and
# FFLAGS. We also set the makefile variable fortran_module_path to a function
# that can be use in the makefile to set where modules are placed or
# searched for.

AC_DEFUN([AC_COMPILER_FLAG],[

#----------------------------------------------------------------
# Set up compiler flags for each supported compiler, handling debug
# and optimized versions differently.
#----------------------------------------------------------------

AC_MSG_CHECKING([whether to enable fortran debug flags])
AC_ARG_ENABLE(debug,
AS_HELP_STRING([--enable-debug],[Enable compiler debugging flags]),
[if test "$enableval" = yes; then  
   enable_debug="yes"
   AC_MSG_RESULT([yes])
   CXXFLAGS="-ggdb -Wall -DBZ_DEBUG -Wno-unused-local-typedefs"
   CFLAGS="-ggdb -Wall"
   if test "$use_ifort" = "yes"; then
     FCFLAGS="-g -O0 -check all,noarg_temp_created -traceback -Difort -heap-arrays 1024"
     FFLAGS="-g -O0  -check all,noarg_temp_created -traceback -Difort -heap-arrays 1024"
   elif test "$use_absoft" = "yes"; then
     FCFLAGS="-g -ggdb -trap=INVALID -trap=DIVBYZERO -trap=OVERFLOW -trap=UNDERFLOW"
     FFLAGS="-g -ggdb -trap=INVALID -trap=DIVBYZERO -trap=OVERFLOW -trap=UNDERFLOW"
   elif test "$use_nag" = "yes"; then
     FCFLAGS="-g -C=all -gline -mtrace=all"
     FFLAGS="-g -C=all -gline -mtrace=all"
   elif test "$use_g95" = "yes"; then
     FCFLAGS="-g -fbounds-check -ftrace=full -freal=nan -fpointer=invalid"
     FFLAGS="-g -fbounds-check -ftrace=full -freal=nan -fpointer=invalid"
   elif test "$use_gfortran" = "yes"; then
     FCFLAGS="-g -fbounds-check -Wall -Wextra"
     FFLAGS="-g -fbounds-check -Wall -Wextra"
   elif test "$use_pathscale" = "yes"; then
     FCFLAGS="-g -ffortran-bounds-check -Wall"
     FFLAGS="-g -ffortran-bounds-check -Wall"
   elif test "$use_portland" = "yes"; then
     FCFLAGS="-g -Mbounds"
     FFLAGS="-g -Mbounds"
   else
     FCFLAGS="-g"
     FFLAGS="-g"
   fi
fi],
[AC_MSG_RESULT([no])
   enable_debug="no"
   if test "$use_ifort" = "yes"; then
# See ticket #526 in trac for why we need the -fp-mode options.
     FCFLAGS="$FCFLAGS -xSSE2 -O2 -fp-model precise -fp-model source -Difort -heap-arrays 1024"
     FFLAGS="$FFLAGS -xSSE2 -O2 -fp-model precise -fp-model source -Difort -heap-arrays 1024"
   elif test "$use_absoft" = "yes"; then
     FCFLAGS="$FCFLAGS -O2"
     FFLAGS="$FFLAGS -O2"
   elif test "$use_nag" = "yes"; then
     FCFLAGS="$FCFLAGS -O4"
     FFLAGS="$FFLAGS -O4"
   elif test "$use_g95" = "yes"; then
     FCFLAGS="$FCFLAGS -O2"
     FFLAGS="$FFLAGS -O2"
   elif test "$use_gfortran" = "yes"; then
     FCFLAGS="$FCFLAGS -O2"
     FFLAGS="$FFLAGS -O2"
   elif test "$use_pathscale" = "yes"; then
     FCFLAGS="$FCFLAGS -Ofast"
     FFLAGS="$FFLAGS -Ofast"
   elif test "$use_portland" = "yes"; then
     FCFLAGS="$FCFLAGS -fastsse -Mipa=fast"
     FFLAGS="$FFLAGS -fastsse -Mipa=fast"
   else
     FCFLAGS="$FCFLAGS -O"
     FFLAGS="$FFLAGS -O"
   fi
])

#----------------------------------------------------------------
# Add support for C++11
#----------------------------------------------------------------

AX_CXX_COMPILE_STDCXX_11([],[mandatory])

#----------------------------------------------------------------
# Add extra checking for fortran, if requested.
#----------------------------------------------------------------

AC_MSG_CHECKING([whether to enable additional checking for fortran])
AC_ARG_ENABLE(addcheck,
AS_HELP_STRING([--enable-addcheck],[Enable additional checking for compiler. Right now, this just adds checking if we are using ifort]),
[if test "$enableval" = yes; then  
   enable_addcheck="yes"
   AC_MSG_RESULT([yes])
   if test "$use_ifort" = "yes"; then
     FCFLAGS="$FCFLAGS -check all,noarg_temp_created"
     FFLAGS="$FFLAGS -check all,noarg_temp_created"
   fi
fi],
[AC_MSG_RESULT([no])
   enable_addcheck="no"
])

#----------------------------------------------------------------
# Determine module information, and change a few things for ifort
#----------------------------------------------------------------

if test "$use_ifort" = "yes"; then
# ifort wants to use its own archive tool for static libraries
   AC_SUBST([AR], [xiar])
# The -Xcompiler is a libtool flag. It tells it not to interpret the -module
# flag as a directive to create a shared modules (e.g. a file that can be 
# dlopened), but rather  pass it on to the compiler without change (which 
# uses is for a Fortran module). Unfortunate clash between two very different
# uses of the term "module" that happen to use the same flag.
   AC_SUBST([fortran_module_path], ['-Xcompiler -module -Xcompiler $(1)'])
# Add -R to options if set in LIBRARY_PATH set by ifort setup

if test "${LIBRARY_PATH}" != ""; then
   AC_SUBST([fortran_extra_ldflags], ["-R ${LIBRARY_PATH}"])
fi
elif test "$use_absoft" = "yes"; then
   AC_SUBST([fortran_module_path], ['-p$(1)'])
elif test "$use_nag" = "yes"; then
   AC_SUBST([fortran_module_path], ['-I$(1) -mdir $(1)'])
elif test "$use_g95" = "yes"; then
   AC_SUBST([fortran_module_path], ['-I$(1) -fmod=$(1)'])
elif test "$use_gfortran" = "yes"; then
   AC_SUBST([fortran_module_path], ['-I$(1) -J$(1)'])
elif test "$use_pathscale" = "yes"; then
   AC_SUBST([fortran_module_path], ['-I$(1) -Xcompiler -module -Xcompiler $(1)'])
elif test "$use_portland" = "yes"; then
   AC_SUBST([fortran_module_path], ['-Xcompiler -module -Xcompiler $(1)'])
else
   AC_SUBST([fortran_module_path], [''])
fi

])
