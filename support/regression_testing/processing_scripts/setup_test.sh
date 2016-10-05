#!/bin/sh

if [ -z "$1" ]; then
    echo "Usage $0 <L2_fp executable>"
    exit 1
fi

# Die on error
set -e

rundir=`pwd`
rundir=`readlink -f $rundir`

basedir=`cd \$( dirname $0 ); pwd`
basedir=`readlink -f $basedir`

export L2_DATASETS_PATH=/groups/algorithm/l2_fp/datasets

. $basedir/common.sh
for dstdir in "${run_list[@]}"; do
    if [ ! -d $dstdir ]; then
	srcdir=$basedir/$dstdir
	mkdir -p $dstdir
	cp -av $srcdir/*.config $srcdir/source_files $dstdir
    fi

    cd $dstdir
    populate.py -b $1 *.config
    cd ..
done

if [ "$rundir" != "$basedir" ]; then
    cp $basedir/common.sh .
    cp $basedir/launch_jobs.sh .
    cp $basedir/run_results.sh .
    cp $basedir/compare_regression_test.sh .
    cp $basedir/splice_results.sh .

    sed -i "s|testdir=.*|testdir=$basedir|" compare_regression_test.sh

    chmod +x launch_jobs.sh run_results.sh compare_regression_test.sh
fi

echo "Need to run on fullerene. To run, you'll do:"
echo "  ./launch_jobs.sh -q long"
echo ""
echo "Check results with ./run_results.sh, timing_summary, and "
echo "./compare_regression_test.sh"
