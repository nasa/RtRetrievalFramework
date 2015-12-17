#!/usr/bin/env bash

script_dir=`cd \$( dirname $0 ); pwd`
run_dir=$PWD
source $script_dir/../inputs_config.sh

if (( $# < 2 )); then
    echo "Not enough arguments"
    exit 1
fi

dataset_name=$1
id_list=$2

if [ ! -e $id_list ]; then
    echo "sounding id file: $id_list does not exist"
    exit 1
fi

base_output_dir=`pwd`
eval dataset_output_dir=$OUTPUT_FILES_STRUCTURE

echo "Assembling dataset $dataset_name at $dataset_output_dir"

mkdir -p $dataset_output_dir
cd $dataset_output_dir

# Extract days from sounding ids to use in find command
days_tmp=`mktemp`
cat $id_list | sed -r 's/^[0-9]{2}//' | sed -r 's/[0-9]{6}$//g' | sort | uniq > $days_tmp

echo "-- Finding $dataset_name input files"

# Create a file to track missing data
missing_data="missing.txt"
echo -n > $missing_data

splice_files_tmp1=`mktemp`
echo > $splice_files_tmp1
for day in `cat $days_tmp`; do
    echo "Searching for $day files"
    orbits=`find /acos/product/Production/v${GOSAT_DATA_VERSION}/L1b${SDOS_PROCESSING_VERSION}/*/${day}/ -name "*.h5" 2>/dev/null | sed -r 's/.*_([0-9]{2})_Production.*/\1/'`
    for curr_orbit in $orbits; do
        echo "Orbit: $curr_orbit"
        missing_idx=0
        curr_splice_files=""
        for file_type in L1b Ecm ; do
            found_file=`find /acos/product/Production/v${GOSAT_DATA_VERSION}/${file_type}${SDOS_PROCESSING_VERSION}/*/${day} -name "acos_${file_type}_${day}_${curr_orbit}_*.h5" 2>/dev/null | head -n 1`
            if [ ! -e "$found_file" ]; then
                # Print to stderr missing file types in such a way that they
                # all are printed on the same line
                if [ "$missing_idx" -eq 0 ]; then
                    echo >> $missing_data
                    echo -n "$day missing " >> $missing_data 
                fi
                echo -n "$file_type " >> $missing_data
                missing_idx=`expr $missing_idx + 1`
            else
                curr_splice_files="$curr_splice_files $found_file"
            fi
        done
        file_count=`echo $curr_splice_files | tr ' ' '\n' | wc -l`
        if [ "$file_count" -lt 2 ]; then
            echo "Missing all necessary files for $day orbit $curr_orbit, only found $file_count file(s)"
            echo "->>$curr_splice_files<<-"
        else
            echo $curr_splice_files >> $splice_files_tmp1 
        fi
    done
done
rm $days_tmp

# Remove any empty spaces
splice_files_tmp2="${run_dir}/${dataset_output_dir}/input_files.txt"
grep -v "^$" $splice_files_tmp1 | sort | uniq > $splice_files_tmp2
rm $splice_files_tmp1

input_list_count=`echo $splice_files_tmp2 | tr ' ' '\n' | wc -l`
if [ "$input_list_count" -gt 0 ]; then
    splice_acos_hdf_files.py -i $splice_files_tmp2 -s $id_list -o acos_L1b${SDOS_PROCESSING_VERSION}_${dataset_name}.h5 -a acos_Ecm${SDOS_PROCESSING_VERSION}_${dataset_name}.h5 

    cd $base_output_dir

    echo "-- Creating config file for inputs"
    . $script_dir/../common_scripts/create_dataset_config.sh
else
    echo "Input files list empty for $dataset_name"
fi
