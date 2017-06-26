#!/bin/sh

run_srch='oco_l2.run'
max_depth=4

for dir in $*
do
    if [ -d $dir ]
    then
	dir=`readlink -f $dir`
	bname=`basename $dir`

	if [[ "$bname" == "." ]]
	then
	    prefix=""
	else
	    prefix="${bname}_"
	fi

	all_list="${prefix}all.list"
	run_list="${prefix}run.list"
	results_txt="${prefix}results.txt"

	if [ -e $all_list ]
	then
	    echo "Checking $bname"
	    echo "Using existing $all_list"

	    run_results.py -r $all_list -d never_ran -o $run_list -s $results_txt
	    echo "Created:" $run_list
	    cat $results_txt
	else
	    echo "Searching $bname"
	    dir=`echo $dir | sed 's|^\./||'`
	    find $dir -maxdepth $max_depth -name $run_srch -exec dirname {} \; > $all_list
	    echo "Created:" $all_list
	fi	
    else
	echo "$dir is not a directory"
    fi
done
