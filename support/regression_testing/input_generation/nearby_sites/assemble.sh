#/usr/bin/env bash

# Exit on first error to occur
set -e

script_dir=`cd \$( dirname $0 ); pwd`
source $script_dir/../inputs_config.sh

base_output_dir=`pwd`
for site_config in `ls $script_dir/site_configs/*.sh`; do
    dataset_name=`basename $site_config | sed 's|\.sh$||'`
    eval dataset_output_dir=$OUTPUT_FILES_STRUCTURE

    echo "Assembling dataset $dataset_name at $dataset_output_dir"

    mkdir -p $dataset_output_dir
    cd $dataset_output_dir

    echo "-- Finding nearby L1B files"
    $script_dir/scripts/find_nearby_l1b.sh $site_config

    echo "-- Finding soundings closest to site"
    $script_dir/scripts/find_closest_to_site.sh $site_config

    echo "-- Splicing site HDF files"
    $script_dir/scripts/splice_hdf.sh $site_config

    cd $base_output_dir

    echo "-- Creating config file for inputs"
    . $script_dir/../common_scripts/create_dataset_config.sh
done
