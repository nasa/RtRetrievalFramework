#!/bin/sh
script_dir=`dirname $0`

# Given a file containing a list of L1B filenames will reformat those filenames and 
# display a list of matching ecmwf files for each L1B file

for l1b_file in `cat $*`; do
    l1b_bn=`basename $l1b_file | sed 's|\.hdf||' | sed 's|\.h5||'`
    ecmwf_glob=`echo $l1b_file | sed 's|/L1b[^/]*/|/Ecm*/|' | sed 's|/L1B[^/]*/|/ECM*/|' | sed 's|/l1b/|/ecmwf/|' | sed 's|L1Bv[^_]*_|ECMv*_|' | sed 's|acos_...|acos_*|' | sed 's|_L1b|_*|' | sed -r 's|[0-9]{12}\.h.*$|*.*|'`
    ecmwf_file=`ls $ecmwf_glob | grep -v '.met'`
    cloud_glob=`echo $l1b_file | sed 's|/L1b[^/]*/|/Cld*/|' | sed 's|acos_...|acos_*|' | sed 's|_L1b|_*|' | sed -r 's|[0-9]{12}\.h.*$|*.*|'`
    cloud_file=`ls $cloud_glob | grep -v '.txt'`
    if [ ! -z "$cloud_file" ];
    then
	echo -e "$l1b_file\t$ecmwf_file\t$cloud_file"
    fi
done
