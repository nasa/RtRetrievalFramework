#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage $0 site_config.sh"
    exit 1
fi

# Site configuraton file
source $1

script_dir=`cd \$( dirname $0 ); pwd`
source $script_dir/common.sh

echo "Finding nearest soundings within $search_distance to ${latitude}/${longitude} from ${find_results_file}"
echo "Creating ${closest_results_file}"
closest_soundings.py --lat=${latitude} --lon=${longitude} -d ${search_distance} `cat ${find_results_file}` > ${closest_results_file} 

if [ -z "$DEBUG" ]; then
    rm $find_results_file
fi
