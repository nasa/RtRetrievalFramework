#!/bin/sh
basedir=`dirname $0`
testdir=$basedir
. $basedir/common.sh
for i in "${run_list[@]}"; do
    echo "****** Checking $i ********"
    compare_regression_test $basedir/$i/output $testdir/$i/expected
done
