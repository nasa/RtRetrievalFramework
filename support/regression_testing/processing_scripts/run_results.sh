#!/bin/sh

basedir=`dirname $0`
. $basedir/common.sh
cd $basedir
for i in "${run_list[@]}"; do
    echo "****** Checking $i ********"
    cd $i
    run_results.py
    cd ..
done
