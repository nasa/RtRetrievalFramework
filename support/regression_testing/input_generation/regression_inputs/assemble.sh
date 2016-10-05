#!/usr/bin/env bash

script_dir=`cd \$( dirname $0 ); pwd`

for id_file in $script_dir/sounding_ids_*.list; do
    dataset_name=$(echo $(basename $id_file) | sed 's/sounding_ids_//' | sed 's/.list//')
    $script_dir/assemble_dataset.sh $dataset_name $id_file
done
