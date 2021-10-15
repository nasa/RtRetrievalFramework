#!/usr/bin/env bash
#PBS -l select=ncpus=16

export imap_file="/fake_path/imap.h5"

export met_file="/fake_path/ecmwf.h5"

export spectrum_file="/fake_path/spectrum.h5"

sounding_id_list_filename=create_run_scripts_test/sounding_id.list
l2_agg_fn=create_run_scripts_test/l2_aggregate.h5
l2_plus_more_agg_fn=create_run_scripts_test/l2_plus_more_aggregate.h5
group_size=1
email_address=

export PYTHONPATH=/fake_python_path
export PATH=/fake_bin_path
export LD_LIBRARY_PATH=/fake_lib_path
export LUA_PATH="/fake_path/input/gosat/config/?.lua;/l2_lua_fake_path"

function email_results {
# Right now, the second results email often isn't sent. I don't know why that
# happens, there is no error message or anything. But for now we'll send a
# separate email saying we are done, and then send results separate.
    echo "Email that jobs is finished"
    cat <<EOF | mail -s "Aggregation finished" $email_address
Aggregation has finished for the file
${l2_agg_fn}
EOF
    echo "Results:"
    cat <<EOF
Aggregation has finished for the file
${l2_agg_fn}

File information:
$(h5ls ${l2_agg_fn}/RetrievalHeader/sounding_id_reference 2>&1)

Plus for file information:
$(h5ls ${l2_plus_more_agg_fn}/RetrievalHeader/sounding_id_reference 2>&1)
EOF
    echo "Emailing results"
    cat <<EOF | mail -s "Aggregation results" $email_address
Aggregation has finished for the file
${l2_agg_fn}

File information:
$(h5ls ${l2_agg_fn}/RetrievalHeader/sounding_id_reference 2>&1)

Plus for file information:
$(h5ls ${l2_plus_more_agg_fn}/RetrievalHeader/sounding_id_reference 2>&1)
EOF
    echo "Done with email for aggregation results"
}

if [ ! -z $email_address ]; then
    trap email_results EXIT
fi

if [ -w "/state/partition1/" ]; then
    worker_temp=/state/partition1/
else
    worker_temp=create_run_scripts_test
fi

if [ ! -z $email_address ]; then
   cat <<EOF | mail -s "Aggregation started" $email_address
Aggregation has started for the file
${l2_agg_fn}
EOF
fi

use_subdir=False
do_aggregate=False

# Aggregate all single sounding output hdf files into a single hdf file
if [ ! -e "$l2_agg_fn" ]; then
    # Use find instead of a glob because there could be too much files that
    # could exhaust the command line length limit
    input_files_tmp=$(mktemp)
    if [ "$do_aggregate" = "True" ]; then
       find create_run_scripts_test/output -name "l2_aggregate_*.h5" > $input_files_tmp
    else
       find create_run_scripts_test/output -name "*.h5" | grep -v "l2_aggregate" | grep -v "l1_aggregate" > $input_files_tmp
    fi

    # Make a version of the sounding id file with newlines instead of spaces
    # for use by splice tool
    in_snd_id_tmp=$(mktemp)
    cat $sounding_id_list_filename | tr ' ' '\n' > $in_snd_id_tmp

    echo "Aggregating L2 output files"
    log_file=agg_$(basename ${l2_agg_fn} | sed 's/\.h5$/.log/')
    /l2_support_fake_path/utils/splice_product_files.py --single-file-type -o $l2_agg_fn -i $input_files_tmp -s $in_snd_id_tmp -l create_run_scripts_test/log/$log_file $* -w 16 --temp $worker_temp

    rm $input_files_tmp
    rm $in_snd_id_tmp
else
    echo "L2 aggregated file exists, skipping creation"
fi

if [ ! -e "$l2_plus_more_agg_fn" ]; then
    echo "Creating L2 plus more aggregated file"
    # Extract sounding ids, and reformat into a format that can be
    # used by splice tool Take sounding ids from L2 aggregate so
    # we know which soundings actually completed Otherwise we will
    # have empty values for L2 datasets when we include IMAP or
    # ABO2 data.

    l2_snd_id_tmp1=$(mktemp)
    l2_snd_id_tmp2=$(mktemp)
    h5dump --noindex -o $l2_snd_id_tmp1 -d RetrievalHeader/sounding_id_reference $l2_agg_fn > /dev/null

    # This line below, turns commas into new lines, removes empty
    # spaces and blank lines
    cat $l2_snd_id_tmp1 | sed 's|,|\n|g' | sed -r 's|^[ ]*||' | grep -v '^$' > $l2_snd_id_tmp2
    rm $l2_snd_id_tmp1

    # If we already have the l1_aggregate_*.h5 for each of the
    # groups l2_fp_job.sh worked on, just combine those.
    if [ "$do_aggregate" = "True" ]; then
    inp_files_tmp=$(mktemp)
    echo $l2_agg_fn > $inp_files_tmp
    find create_run_scripts_test/output -name "l1_aggregate_*.h5" >> $inp_files_tmp

    log_file=agg_$(basename ${l2_plus_more_agg_fn} | sed 's/\.h5$/.log/')
    /l2_support_fake_path/utils/splice_product_files.py --multiple-file-type --splice-all --rename-mapping --agg-names-filter -o $l2_plus_more_agg_fn -i $inp_files_tmp -s $l2_snd_id_tmp2 -l create_run_scripts_test/log/$log_file $* -w 16 --temp $worker_temp

    rm $inp_files_tmp
    rm $l2_snd_id_tmp2
    else
    # Otherwise, generate everything from scratch
    # Combine L1B and L2 files into one file with IMAP and ABand
    # files if they are supplied
    inp_files_tmp=$(mktemp)
    echo $l2_agg_fn > $inp_files_tmp
    if [ ! -z "$input_file_mapping" ] && [ -e "$input_file_mapping" ]; then
        # Set up input files from file mapping
        while read -r sounding_id file_map; do
            eval $(echo $file_map | tr ';' '\n')
            if [ ! -z "$spectrum_file" ] && [ ! -z "$imap_file" ] && [ ! -z "$aband_file" ]; then
                echo $spectrum_file
                echo $imap_file 
                echo $aband_file 
            elif [ ! -z "$spectrum_file" ] && [ ! -z "$imap_file" ]; then
                echo $spectrum_file
                echo $imap_file 
            fi
        done < $input_file_mapping | sort | uniq >> $inp_files_tmp
    else
        # Use input files from script variables
        for fn in $spectrum_file $imap_file $aband_file; do
            echo $fn >> $inp_files_tmp
        done
    fi

    log_file=agg_$(basename ${l2_plus_more_agg_fn} | sed 's/\.h5$/.log/')
    /l2_support_fake_path/utils/splice_product_files.py --multiple-file-types --splice-all --rename-mapping --agg-names-filter -o $l2_plus_more_agg_fn -i $inp_files_tmp -s $l2_snd_id_tmp2 -l create_run_scripts_test/log/$log_file
    rm $l2_snd_id_tmp2 $inp_files_tmp
    fi
    # Create retrieval_index dataset based on L1B file
    if [ ! -z "$spectrum_file" ]; then
        echo "Adding retrieval information datasets"
        /l2_support_fake_path/utils/add_spliced_retrieval_info.py $spectrum_file $l2_plus_more_agg_fn
    fi
else
    echo "L2 plus more file exists, skipping creation"
fi

