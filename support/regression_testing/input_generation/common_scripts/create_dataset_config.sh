# This should be sourced into the place where the script defines 
# $dataset_name and wants to create the config
#
# Also this expects to be run from the base directory we started at

base_dir=`pwd`

eval config_dir=$CONFIG_FILE_LOCATION
config_dir=`readlink -f $config_dir`

eval data_dir=$OUTPUT_FILES_STRUCTURE
data_dir=`readlink -f $data_dir`

mkdir -p ${config_dir}
cd ${config_dir}

config_filename=${config_dir}/${dataset_name}.config
create_config.py -t gosat ${data_dir}/* -o ${config_filename}

# Make paths in config relative
sed -i "s|${config_dir}/||" ${config_filename}

cd $base_dir
