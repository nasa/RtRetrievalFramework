# Specifies the SDOS processing and GOSAT data versions to use when building 
# regression testing inputs

GOSAT_DATA_VERSION="161160"
SDOS_PROCESSING_VERSION="30502"

# $dataset_name here should be evaluated by the
# assemble.sh script after setting the value
# of $dataset_name
OUTPUT_FILES_STRUCTURE='${dataset_name}/source_files'
CONFIG_FILE_LOCATION='${dataset_name}'
