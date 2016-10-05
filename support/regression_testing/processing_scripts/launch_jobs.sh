#!/bin/sh

if [ -z "$1" ]; then
    echo "Usage $0 -q long"
    exit 1
fi
basedir=`dirname $0`
. $basedir/common.sh
for i in "${run_list[@]}"; do
    $basedir/$i/launch_jobs.sh $*
done
