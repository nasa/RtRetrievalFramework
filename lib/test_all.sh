#!/bin/sh
# You can edit this on your local copy to run a single test or a selected set
# of tests. This can be useful when you are working in one area and want a
# fast turn around. You should however always have all the tests on
# when you check this in.

echo "==========================================================="
echo "Can specify run_test=<exp> on command line to run subset of"
echo "unit tests. See boost unit test documentation for details."
echo ""
echo "If log_test is defined, add verbose information as each"
echo "test is run."
echo ""
echo "If valgrind is defined, then we run the test with valgrind."
echo ""
echo "Do long_check to run all tests, check to run shorter tests"
echo "==========================================================="

if [ ${valgrind} ] ; then
    tool_command="valgrind --max-stackframe=5000000 --error-exitcode=1 --track-origins=yes --suppressions=$(dirname $0)/../config/valgrind.suppressions"
elif [ ${gdb} ]; then
    tool_command="gdb --args"
else
    tool_command=""
fi

if [ ${log_test} ] ; then
    ${tool_command} ./lib/test_all --log_level=test_suite --run_test=${run_test}
else
    ${tool_command} ./lib/test_all --show_progress --run_test=${run_test}
fi

# Valgrind version of test run. 
# Note valgrind *hates* ifort, you'll get lots of fortran errors if you use
# that. If you are going to use valgrind, you should used gfortran-4.5 rather
# than ifort.
# 
# Note thate the version of valgrind in /groups/algorithm/tools/install/bin
# is newer than the one found on the system. You can add  
# --track-origins=yes to track down where things like uninitialized data
# comes from.
#/groups/algorithm/tools/install/bin/valgrind --suppressions=${L2_EXE_PATH}/unit_test_data/valgrind.suppressions ./lib/test_all --log_level=test_suite --run_test=${run_test}
