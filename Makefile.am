## Process this file with automake to produce Makefile.in

#=================================================================
# Based on Peter Millers well known paper "Recursive Make Considered Harmful"
# (at http://miller.emu.id.au/pmiller/books/rmch/), we have the whole system
# described and built with this one Makefile.

# Tell automake to put the object file for foo/bar.c in dir foo
AUTOMAKE_OPTIONS = subdir-objects

# For right now, continue to use serial-tests. We'll want to change this
# in the future
AUTOMAKE_OPTIONS+= serial-tests

# Don't worry about Makefile portability warnings. We can assume that we 
# are using GNU make. Pretty much any modern system will have GNU make on it.

AUTOMAKE_OPTIONS += --warnings=no-portability 

# Local autoconf packages
ACLOCAL_AMFLAGS = -I config/m4

# Don't wrap executable in a libtool script. This is useful if we install
# somewhere other than where we build. I'm not sure about this, we may
# want to actual do the install. But for now, turn this off
#AM_LDFLAGS = -no-install

# Build these up as we go
lib_LTLIBRARIES =
noinst_LTLIBRARIES =
CLEANFILES =
DISTCLEANFILES =
bin_PROGRAMS =
bin_SCRIPTS =
EXTRA_DIST = config/depf90
EXTRA_LTLIBRARIES =
BUILT_SOURCES = 
fullphysicsincdir = @includedir@/full_physics
ourboostincdir = @includedir@/boost
fullphysicsinc_HEADERS = 
ourboostinc_HEADERS = 
lib_test_all_SOURCES =
config_DATA =
input_DATA =
AM_CPPFLAGS =
AM_LDFLAGS = 
phony =
INSTALL_DATA_HOOK =
pkgpython_LTLIBRARIES = 
pkgpython_PYTHON =
SWIG_SRC =
SWIG_FLAG =

# Make sure to look in local install directory first, in case we are
# building our own copy of libraries also found in thirdparty library area.
# However boost/backtrace is a special case, we want this before include
# otherwise we get a build loop
AM_CPPFLAGS+= -I$(srcdir)/@supportsrc@
AM_CPPFLAGS+= -I$(prefix)/include
AM_LDFLAGS+= -L$(prefix)/lib
AM_LDFLAGS+= $(fortran_extra_ldflags)

# Include boost in builds. We have two flavors, one that uses shared libraries
# and one that uses static libraries.
AM_CPPFLAGS+= $(BOOST_CPPFLAGS)
BOOST_LDFLAGS = -L$(BOOST_LIBDIR) -R $(BOOST_LIBDIR) 
BOOST_LDFLAGS+= -lboost_regex -lboost_date_time -lboost_iostreams 
BOOST_LDFLAGS+= -lboost_filesystem -lboost_system
if HAVE_BOOST_STATIC
  BOOST_LDSTATIC = $(BOOST_LIBDIR)/libboost_regex.a
  BOOST_LDSTATIC+= $(BOOST_LIBDIR)/libboost_date_time.a
  BOOST_LDSTATIC+= $(BOOST_LIBDIR)/libboost_iostreams.a
  BOOST_LDSTATIC+= $(BOOST_LIBDIR)/libboost_filesystem.a
  BOOST_LDSTATIC+= $(BOOST_LIBDIR)/libboost_system.a
endif

# Include blitz in builds
AM_CPPFLAGS+= $(BLITZ_CPPFLAGS)

# Include CPPAD in builds (may go away)
AM_CPPFLAGS+= $(CPPAD_CPPFLAGS)

# Include GSL in builds
AM_CPPFLAGS+= $(GSL_CPPFLAGS)

# We don't actually use any threading code, but on some systems the static
# libraries such as libgfortran require pthread library. If we found these
# libraries during configuration, then go ahead and include on the link line.
# If this isn't actually needed, it won't hurt to include it anyways.
if HAVE_PTHREAD
AM_LDFLAGS+= $(PTHREAD_LIBS)
endif

#=================================================================
# Include directions for building doxgen stuff.

include config/am/doxygen.am
if WITH_DOCUMENTATION
doc/FullPhysics.tag: config/am/doxygen.am config/doxygen.cfg 
doc/FullPhysics.tag: config/doc_footer.html
doc/FullPhysics.tag: $(libradiant_la_SOURCES) $(libfp_la_SOURCES)
doc/FullPhysics.tag: $(l2_fp_SOURCES) $(libfp_la_SOURCES)
doc/FullPhysics.tag: $(FLIB_FILES)
endif

#=================================================================
# Unit test
check_PROGRAMS = lib/test_all
lib_test_all_SOURCES += lib/test_all.cc
lib_test_all_LDADD = libfp.la
# The Mac doesn't support no-as-needed, so we leave it off there.
if IS_MAC
lib_test_all_LDFLAGS = -no-install
else
lib_test_all_LDFLAGS = -no-install -Wl,--no-as-needed
endif
lib_test_all_LDFLAGS+= $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
lib_test_all_LDFLAGS+= $(FCLIBS) $(PTHREAD_LIBS)
TESTS = lib/test_all.sh
EXTRA_DIST+= lib/test_all.sh
EXTRA_DIST+= config/valgrind.suppressions
# Variables used in testing
export abs_top_srcdir
export abscodir
export merradir
export run_test = *
export LUA_PATH = $(abs_srcdir)/input/common/config/?.lua;$(abs_srcdir)/input/gosat/config/?.lua;$(abs_srcdir)/input/oco/config/?.lua;$(abs_srcdir)/input/fts/config/?.lua;$(abs_srcdir)/unit_test_data/?.lua;

#=================================================================
# Normally check builds all the libraries and executables. This
# version builds just what is needed for unit tests, which is faster
# to do during development.

fast_check: $(BUILT_SOURCES)
	$(MAKE) $(AM_MAKEFLAGS) fast_check-am

fast_long_check: $(BUILT_SOURCES)
	$(MAKE) $(AM_MAKEFLAGS) L2_FP_LONG_CHECK=1 fast_check-am

fast_timing_check: $(BUILT_SOURCES)
	$(MAKE) $(AM_MAKEFLAGS) L2_FP_TIMING_CHECK=1 fast_check-am

fast_check-am: lib/test_all
	@echo "====================================================="
	@echo "This just builds and runs the unit tests (useful during"
	@echo "development). It does not rebuild the executables - use"
	@echo "'make check' for that."
	@echo "====================================================="
	$(srcdir)/lib/test_all.sh

phony += check_message long_check fast_long_check fast_timing_check
phony += fast_check timing_check nosetests_check

.PHONY: $(phony)

check: check_message

check: nosetests_check

check_message:
	@if test "${L2_FP_LONG_CHECK}" = ""; then \
          echo "====================================================="; \
          echo "Note that there are longer running tests that aren't"; \
          echo "run with make check. To run them, do make long_check"; \
          echo "Timing checks are run byu make timing_check";          \
          echo ""; \
          echo "You can also run fast_check, which only builds the"; \
          echo "unit tests, not executables (useful for development)"; \
          echo "====================================================="; \
        fi

long_check:
	@echo "====================================================="
	@echo "This does extra tests which take a longer time to run"
	@echo "(in addition to all the short tests)."
	@echo "====================================================="
	$(MAKE) L2_FP_LONG_CHECK=1 check

# Temporary, turn this off because it fails in buildbot. We'll add this
# back in when we have whatever issue worked out.
#	$(MAKE) L2_FP_LONG_CHECK=1 check && \
#        $(MAKE) pythondoctest

really_long_check:
	@echo "====================================================="
	@echo "This does extra tests which take a longer time to run"
	@echo "(in addition to all the short tests)."
	@echo "This runs stuff that takes a significant amount of time"
	@echo "to run and that we don't normally need to do. You"
	@echo "don't normally need to run these tests"
	@echo "====================================================="
	$(MAKE) L2_FP_REALLY_LONG_CHECK=1 check

timing_check:
	@echo "====================================================="
	@echo "This does timing tests (in addition to short tests)"
	@echo "It can be useful to run this with run_test=*/*timing"
	@echo "to just run timing tests"
	@echo "====================================================="
	$(MAKE) L2_FP_TIMING_CHECK=1 check

#=================================================================
# Run nosetests

if HAVE_NOSETESTS
nosetests_check: install-python
	export PYTHONPATH=$(pythondir):$(pyexecdir):$(PYTHONPATH) && \
        $(NOSETESTS) $(abs_srcdir)/lib/Python $(abs_srcdir)/verification_test

else # HAVE_NOSETESTS
nosetests_check:
	@echo "====================================================="
	@echo "Skipping nosetests because it wasn't found on the"
	@echo "system during configuration"
	@echo "====================================================="
endif

#=================================================================
# Include source files
include lib/lib.am
include bin/bin.am
include thirdparty/thirdparty.am
include script/script.am
include support/support.am

#=================================================================
# Configuration input for full physics

include input/input.am

#=================================================================
# Bindings for python from Swig

include bindings/python/python.am

#=================================================================
# Copy diagrams for documentation

doxygen_run_extra += doc/html/overall.png

doc/html/overall.png: design/overall.png
	cp $^ $@

#=================================================================
# Extra files to include in the distribution
EXTRA_DIST += nonjpl_build.script README.developer
EXTRA_DIST += setup_env.sh
EXTRA_DIST += config/setup_compiler_env.sh
EXTRA_DIST += config/doc_footer.html config/doxygen.cfg 
EXTRA_DIST += config/doxygen_fortran.cfg
EXTRA_DIST += doc/LIDORT_userguide_v3.0.pdf
EXTRA_DIST += doc/users_guide.pdf
EXTRA_DIST += input/setup_env.sh
EXTRA_DIST += input/gosat
EXTRA_DIST += operations support unit_test_data

#=================================================================
# Test at a higher level than unit tests

include test/test.am

#=================================================================
# Sphinx python documentation

include doc/python/doc_python.am

install-data-hook: $(INSTALL_DATA_HOOK)

