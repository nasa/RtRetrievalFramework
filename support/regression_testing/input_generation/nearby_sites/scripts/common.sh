# Common varialbes used by the scripts here

script_dir=`cd \$( dirname $0 ); pwd`

# This is where GOSAT_DATA_VERSION and SDOS_PROCESSING_VERSION
# are defined
source $script_dir/../../inputs_config.sh

# These are temporary files needed between the scripts
find_results_file="${site_name}_l1b_v${SDOS_PROCESSING_VERSION}_found.list"
closest_results_file="${site_name}_v${SDOS_PROCESSING_VERSION}_closest.txt"
