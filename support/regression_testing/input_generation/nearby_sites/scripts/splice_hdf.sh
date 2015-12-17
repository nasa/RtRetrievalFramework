#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage $0 site_config.sh"
    exit 1
fi

# Site configuraton file
source $1

script_dir=`cd \$( dirname $0 ); pwd`
source $script_dir/common.sh

closest_l1b_file="${site_name}_l1b_v${SDOS_PROCESSING_VERSION}_closest.list"
closest_addl_file="${site_name}_l1b_plus_addl_v${SDOS_PROCESSING_VERSION}_closest.list"
sounding_id_file="${site_name}_v${SDOS_PROCESSING_VERSION}_sounding_id.list"

grep L1b ${closest_results_file} > ${closest_l1b_file}
grep ^200 ${closest_results_file} | awk '{print $1}' > ${sounding_id_file}

echo "Finding associated ECMWF and Cloud files"
find_l1b_matching_ecmwf_cloud.sh ${closest_l1b_file}  > ${closest_addl_file}

splice_acos_hdf_files.py -i ${closest_addl_file} -s ${sounding_id_file} -o ${site_name}_l1b_${SDOS_PROCESSING_VERSION}_spliced.h5 -a ${site_name}_ecmwf_${SDOS_PROCESSING_VERSION}_spliced.h5 -a ${site_name}_cloud_${SDOS_PROCESSING_VERSION}_spliced.h5 

if [ -z "$DEBUG" ]; then
    rm $closest_results_file $closest_l1b_file $closest_addl_file $sounding_id_file
fi
