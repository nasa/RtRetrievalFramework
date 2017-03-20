#!/bin/bash

ABSCO_DIR="/groups/algorithm/l2_fp/absco"

function message
{
    echo "$*" >&2
}

function usage
{
    if [ ! -z "$1" ]; then
        message $1

    fi
    message "Usage: $0 <L1B or log filename> <input config file> [<arguments>]"
    message "Supply L1B filename only if log is located in the same directory as the L1B"
    message "Optional arguments must go after the required log and config filenames"
    message ""
    message "  -s <filename>  - Alternative static input filename"
    message "  -b <name>      - Only retrieve for the given band/spectrometer name"
    message "  -o <o2_concen> - O2 concentration to use"
    message "  -m <met_File>  - Analyze to see which ECMWF reader to use"
    message "  -c <co2_scale> - A Prior CO2 scaling"
    message "  -a             - Override aerosol detection and use default" 
    exit 1
}

function error
{
    message $1
    exit 1
}

if [ $# -lt 2 ]; then
    usage "Not enough arguments"
fi

# Try and determine log filename from HDF filename if that is what is supplied
log_or_hdf=$1
if [[ $log_or_hdf =~ \.hdf$ || $log_or_hdf =~ \.h5$  ]]; then
    log_name=`basename $log_or_hdf | sed s/\.h.*$/.log/`
    hdf_dirname=`dirname $log_or_hdf`
    hdf_realdir=`readlink -f $log_or_hdf`
    hdf_realdir=`dirname $hdf_realdir`

    if [ -e ${hdf_dirname}/${log_name} ]; then
        log_filename=${hdf_dirname}/${log_name}
    elif [ -e ${hdf_realdir}/${log_name} ]; then
        log_filename=${hdf_realdir}/${log_name}
    else
        error "Could not find log file from $log_or_hdf"
    fi
elif [[ $log_or_hdf =~ \.log$ && -e "$log_or_hdf" ]]; then
    log_filename=$log_or_hdf
else
    error "A valid log filename was not specifcied: $log_or_hdf"
fi

if [ -e "$log_filename" ]; then
    message "Found log filename: $log_filename"
else
    error "Log filename does not exist: $log_filename"
fi

# Make sure we have a good config filename
config_filename=$2
if [ ! -e "$config_filename" ]; then
    error "Config file does not exist: $config_filename"
elif [[ ! $config_filename =~ \.lua$ ]]; then
    error "Config file does not have .lua extension: $config_filename"
fi

# Process optiona arguments
OPTIND=3
while getopts  "b:s:o:m:c:a" flag
do
    if [[ "$flag" == "b" ]]; then
        # Optional specification of a particular band to retrieve
        which_band=$OPTARG
    elif [[ "$flag" == "s" ]]; then
        # Optional specification of a particular static file
        static_file=$OPTARG
    elif [[ "$flag" == "o" ]]; then
        # Optionally change O2 concentration
        o2_concentration=$OPTARG
    elif [[ "$flag" == "m" ]]; then
        met_file=$OPTARG
    elif [[ "$flag" == "c" ]]; then
        ap_co2_scaling=$OPTARG
    elif [[ "$flag" == "a" ]]; then
        dflt_aerosol=1
    fi
done

# Look through log file and build up options
config_options_file=`mktemp`
message "Writing config options to temporary file: $config_options_file"

echo > $config_options_file
echo "-- == Configuration deduced from $log_filename ==" >> $config_options_file

# Configure ABSCO version
absco_version=`grep 'MS3 absorption coefficient file prefix' $log_filename | cut -d "'" -f 2 | sed -r 's/^.*_(v[0-9]+)(_unscaled)?_ms3.*$/\1\2/'` 
absco_ver_dir=""
absco_co2_fn=""
absco_h2o_fn=""
absco_o2_fn=""

if [[ "$absco_version" == "v411_unscaled" ]]; then
    absco_co2_fn="v4.1.1/lowres/co2_v4.1.1-lowres.hdf"
    absco_h2o_fn="v4.1.1/lowres/h2o_v4.1.1-lowres.hdf"
    absco_o2_fn="v4.1.1/lowres/o2_v4.1.1-lowres.hdf"
    co2_scale="1.0"
    o2_scale="1.0"
elif [[ "$absco_version" == "v411" ]]; then
    absco_co2_fn="v4.1.1_rescaled/lowres/co2_v4.1.1-lowres-sco2_rescaled_by_0.99.hdf"
    absco_h2o_fn="v4.1.1_rescaled/lowres/h2o_v4.1.1-lowres.hdf"
    absco_o2_fn="v4.1.1_rescaled/lowres/o2_v4.1.1-lowres-rescaled_by_1.0125.hdf"
    co2_scale="{1.0, 1.0038, 0.9946}"
    o2_scale="1.0125"
elif [[ "$absco_version" == "v420_unscaled" ]]; then
    absco_co2_fn="v4.2.0_unscaled/co2_v4.2.0_with_ctm.hdf"
    absco_h2o_fn="v4.2.0_unscaled/h2o_v4.2.0.hdf"
    absco_o2_fn="v4.2.0_unscaled/o2_v4.2.0_drouin.hdf"
    co2_scale="1.0"
    o2_scale="1.0"
else
    error "Did not know how to use ABSCO version string: $absco_version"
fi

cat >>$config_options_file <<OPTS

-- ABSCO paths and filenames
config.fm.atmosphere.absorber.CO2.absco = "$absco_co2_fn"
config.fm.atmosphere.absorber.H2O.absco = "$absco_h2o_fn"
config.fm.atmosphere.absorber.O2.absco = "$absco_o2_fn"
config.fm.atmosphere.absorber.CO2.table_scale = $co2_scale
config.fm.atmosphere.absorber.O2.table_scale = $o2_scale
OPTS

# Set cloud free RT
# Either aerosols are turned off in MS3

# or the scene file does not include aerosols
if [[ -z "$dflt_aerosol" ]]; then
    cloud_free_opt=`grep 'Cloud-free RT' $log_filename` 
    scene_filename=`grep 'Simulation file name' $log_filename | cut -d "'" -f 2`
    if [[ $cloud_free_opt =~ T$ || $scene_filename =~ T_T_T ]]; then
        cat >>$config_options_file <<OPTS

-- Turn off aerosol retrieval and scattering calculations
config.fm.atmosphere.aerosol.creator = ConfigCommon.rayleigh_only
OPTS
    fi
fi

# Turn off solar doppler shift if off in the simulator
if [[ `grep 'Use solar Doppler shift' $log_filename` =~ F$ ]]; then
    cat >>$config_options_file <<OPTS

-- Turn off solar doppler shift
config.fm.spectrum_effect.solar_model.doppler_shift.do_doppler_shift = false 
OPTS
fi

# Turn off instrument doppler shift if disabled in simulator
if [[ `grep 'Use instrument Doppler shift' $log_filename` =~ F$ ]]; then
    cat >>$config_options_file <<OPTS

-- Turn off instrument doppler by removing this item from the spectrum effect list
config.fm.spectrum_effect.speceff["instrument_doppler"] = nil
OPTS
fi

# Turn on fluorescence if included in orbit simulator
if [[ `grep 'Turn on plant fluorescence' $log_filename` =~ F$ ]]; then
    cat >>$config_options_file <<OPTS

-- Disable fluorescence
config.fm.spectrum_effect.speceff["fluorescence"] = nil
for i,k in ipairs(config.fm.spectrum_effect.speceff) do
    if k == "fluorescence" then
        table.remove(config.fm.spectrum_effect.speceff, i)
    end
end
OPTS
fi

# Write out ILS extent configuration
ils_extents=`grep 'ILS extent' $log_filename | sed -r 's/^# ILS extent\s+//'`
unit_code=`grep Units $log_filename | sed 's/.*\([0-9]\)$/\1/'`
if [ $unit_code == 0 ]; then
    units="cm^-1"
elif [ $unit_code == 1 ]; then
    units="um"
elif [ $unit_code == 2 ]; then
    units="nm"
else
    error "Unknown units code: $unit_code"
fi
cat >>$config_options_file << OPTS

-- Match the same ILS extents as used in the simulator
config.fm.instrument.ils_half_width = { 
OPTS
echo -n "    " >> $config_options_file
for curr_ext in $ils_extents; do
    half_ext=`perl -e "printf \"%.2e\", \"$curr_ext\" / 2.0"`
    echo -n "DoubleWithUnit($half_ext, \"$units\"), " >> $config_options_file
done
echo >> $config_options_file
echo "}" >> $config_options_file

# Optionally output setup for a particular band to be retrieved
if [ ! -z "$which_band" ]; then
    cat >>$config_options_file << OPTS

-- Retrieve only band $which_band
config.which_spectrometers = "$which_band"
require "single_band_support"
init_single_band_support(config)
OPTS
fi

# Optionally use a different static file
if [ ! -z "$static_file" ]; then
    static_file=`readlink -f $static_file`
    cat >>$config_options_file << OPTS

-- Use a different static file
config.static_file = "$static_file"
OPTS
fi

# Optionally use a different O2 concentration 
if [ ! -z "$o2_concentration" ]; then
    cat >>$config_options_file << OPTS

-- Use a different O2 concentration 
function user_o2_concentration()
    local o2_val = Blitz_double_array_1d(1)
    o2_val:set(0, $o2_concentration)
    return o2_val
end
config.fm.atmosphere.absorber.O2.apriori = user_o2_concentration
OPTS
fi

# Add use of orbit simulator met file reader
if [ ! -z "$met_file" ]; then
    if [ ! -z "$(h5ls $met_file | grep '^ecmwf')" ]; then
        cat >>$config_options_file << OPTS

-- Use CSU ECMWF format reader
config.fm.input.metg.creator = OcoConfig.oco_meteorology
OPTS
    fi
fi

# Add rescaling of A Priori CO2
if [ ! -z "$ap_co2_scaling" ]; then
        cat >>$config_options_file << OPTS

--- Replace the apriori function with this new one function 
--- with one that scales the original co2 apriori. 
old_aprior_func = config.fm.atmosphere.absorber.CO2.apriori
function config:scaled_apriori()
   local v = old_aprior_func(self)
   for i=0,v:rows() - 1 do
      v:set(i, v(i) * ${ap_co2_scaling})
   end
   return v
end
config.fm.atmosphere.absorber.CO2.apriori = config.scaled_apriori        
OPTS
fi

# Mark end of our changes
echo >> $config_options_file
echo "-- == END ==" >> $config_options_file
echo >> $config_options_file

# Outputs config file now
do_config_line=`grep -n "do_config" $config_filename | cut -d  ':' -f 1`
do_config_line=`expr $do_config_line '-' 1`
sed  "${do_config_line}r $config_options_file" $config_filename

message "Removing temporary file: $config_options_file"
rm $config_options_file
