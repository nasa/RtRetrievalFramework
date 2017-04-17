#!/bin/bash

# Things that vary from one run to the next
export imap_file="/fake_path/imap.h5"

export met_file="/fake_path/ecmwf.h5"

export spectrum_file="/fake_path/spectrum.h5"

export l2_binary_filename=/fake_path/l2_fp
export l2_config_filename=/fake_path/input/gosat/config/config.lua
export processing_dir=create_run_scripts_test
export log_directory=create_run_scripts_test/log
export sounding_id_list_filename=create_run_scripts_test/sounding_id.list
export output_directory=create_run_scripts_test/output
export group_size=1
export PYTHONPATH=/fake_python_path
export PATH=/fake_bin_path
export LD_LIBRARY_PATH=/fake_lib_path
export LUA_PATH="/fake_path/input/gosat/config/?.lua;/l2_lua_fake_path"
if [[ "/fake_absco_path" != "/groups/algorithm/l2_fp/absco" ]]; then
    # Only do this if not the default location. This is because it
    # interfers with the logic of the search for a local copy of the
    # absco file on torque
    export abscodir=/fake_absco_path
fi
export merradir=/fake_merra_path


# The rest of this should not need to be changed.
cd ${processing_dir}
id_list=( `cat ${sounding_id_list_filename}` )

keep_exit="no"
print_help="no"
force_quiet="no"
while getopts :ehqa:m: opt; do
  case $opt in
    a)
      export abscodir=$OPTARG
      ;;
    m)
      export merradir=$OPTARG
      ;;
    e)
      keep_exit="yes"
      ;;
    q)
      force_quiet="yes"
      ;;
    h | \?)
      print_help="yes"
      ;;
  esac
done
shift $((OPTIND-1))
if [[ ! -z "${PBS_ARRAY_INDEX}" ]]; then
   job_index=${PBS_ARRAY_INDEX}
   do_tee="no"
else
   if [[ $# -lt 1 ]]; then
     print_help="yes"
   else
     job_index=$1
     do_tee="yes"
   fi
fi

if [ "$print_help" = "yes" ]; then
cat <<EOF
Usage: $0 [options] [<job_index>]"

Need to supply job array index, or pass in as PBS_ARRAY_INDEX

   -h  Print help
   -a abscodir 
       Supply absco directory to use, rather than one set when
       script 
   -a merradir 
       Supply merra directory to use, rather than one set when
       script 
   -e  Keep exit code from l2 code
   -q  Don't write output to stdout, even if we aren't in PBS
EOF
  exit 1
fi

if [ "$force_quiet" = "yes" ]; then
   do_tee="no"
fi
export do_tee



# Break up log and output directories so they don't get too big. Just
# use the job index padded with 0's for this. Skip if we aren't grouping
# data into large enough groups - we don't want directories with a single
# file in them.
if [ $group_size -gt 9 ]; then
    log_directory=${log_directory}/$(printf "%04d" $job_index)
    output_directory=${output_directory}/$(printf "%04d" $job_index)
fi
mkdir -p ${log_directory}
mkdir -p ${output_directory}

run_job_index() {
    id_index=$1
    id_list=( `cat ${sounding_id_list_filename}` )
    export sounding_id="${id_list[${id_index}]}"

    if [ ! -z "$input_file_mapping" ] && [ -e "$input_file_mapping" ]; then
        # If there is an input file mapping, then extract declarations of filenames from it
        # Expected format is 1 line per sounding id: 
        # Sounding id with space seperating it from key=value pairs seperated by semicolons, ie:
        # <sounding_id> inp_file1=/path/to/file1;inp_file2=/path/to/file2
        # The keywords names are evaluated as variables so they should match the expected enviromental variable names
        # used by this script.
        file_decls=$(grep ^$sounding_id $input_file_mapping | sed "s/^$sounding_id //" | tr ';' '\n' | sed 's/^/export /g')
        eval $file_decls
    fi
    
    # The mac doesn't have the same options in the time command as on linux
    if test ! -z "`uname | grep Darwin`"
    then
       time_command="/usr/bin/time"
    else
       time_command="/usr/bin/time -v -o ${log_directory}/runtime_${sounding_id}.txt"
    fi
    if test -e ${output_directory}/l2_${sounding_id}.h5; then
       echo "Skipping ${sounding_id} because output file already exists"
    else
       if [ "$do_tee" = "yes" ]; then
           ($time_command ${l2_binary_filename} ${l2_config_filename} ${output_directory}/l2_${sounding_id}.h5) 2>&1 | tee ${log_directory}/l2_${sounding_id}.log.running
           l2_run_status=${PIPESTATUS[0]}
       else
           ($time_command  ${l2_binary_filename} ${l2_config_filename} ${output_directory}/l2_${sounding_id}.h5) > ${log_directory}/l2_${sounding_id}.log.running 2>&1
           # In torque, the fact that this ran to complete is "success", even
           # if sounding failed.
           l2_run_status=$?
       fi
       mv -f ${log_directory}/l2_${sounding_id}.log.running ${log_directory}/l2_${sounding_id}.log
    fi
}
export -f run_job_index

beg_id_index=$(expr ${job_index} '*' ${group_size})
end_id_index=$(expr ${beg_id_index} '+' ${group_size} '-' 1)
if [ ${end_id_index} -gt ${#id_list[@]} ]; then
    end_id_index=$(expr ${#id_list[@]} - 1)
fi

seq ${beg_id_index} ${end_id_index} | parallel --gnu -j 1 --ungroup --env run_job_index 'run_job_index {}'

if [ "$keep_exit" = "yes" ]; then
   exit ${l2_run_status}
fi
