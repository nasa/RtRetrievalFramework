#!/bin/bash
basedir=`dirname $0`
. $basedir/common.sh
for run_dir in "${run_list[@]}"; do
    cd $basedir/$run_dir
    if [ -e expected ]; then
	splice_acos_hdf_files.py --aggregate -o l2_${run_dir}_output.h5 expected/*.h5
    else
	splice_acos_hdf_files.py --aggregate -o l2_${run_dir}_output.h5 output/*.h5
    fi
    cd -
done

