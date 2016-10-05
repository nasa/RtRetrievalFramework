#!/usr/bin/env bash

if [ -z "$1" ]; then
    echo "Usage $0 site_config.sh"
    exit 1
fi

# Site configuraton file
source $1

script_dir=`cd \$( dirname $0 ); pwd`
source $script_dir/common.sh

echo "Creating L1B search results file: ${find_results_file}"
echo -n > ${find_results_file}

echo "Searching GOSAT ver: ${GOSAT_DATA_VERSION}"
for curr_path in $paths; do
    echo "Searching orbit path: ${curr_path}"
    find /acos/product/Production/v${GOSAT_DATA_VERSION}/L1b${SDOS_PROCESSING_VERSION}/ -name "acos_L1b_*_${curr_path}_*_L1b${SDOS_PROCESSING_VERSION}_*" >> ${find_results_file}
done



