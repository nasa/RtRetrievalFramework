#!/bin/sh

# Extracts only the essential datasets from a L1B file for a given set of sounding ids
# This script uses splice_product_files.py to accomplish this task and mainly exists
# just to easily se tup the inputs to achieve this.

function usage
{
    echo "Usage: $0 <input_filename> [<sounding_id_1> <sounding_id_2>...]"
    exit 1
}
[ $# -gt 0 ] || usage

datasets_input_file=`mktemp`

cat > ${datasets_input_file} <<DATASETS
/ABandCloudScreen/cloud_flag
/ecmwf/specific_humidity
/ecmwf/specific_humidity_pressures
/ECMWF/specific_humidity_profile_ecmwf
/ecmwf/surface_pressure
/ECMWF/surface_pressure_ecmwf
/ecmwf/temperature
/ecmwf/temperature_pressures
/ECMWF/temperature_profile_ecmwf
/ECMWF/vector_pressure_levels_ecmwf
/ecmwf/windspeed_u
/ECMWF/windspeed_u_ecmwf
/ecmwf/windspeed_v
/ECMWF/windspeed_v_ecmwf
/Meteorology/ozone_profile_met
/Meteorology/specific_humidity_profile_met
/Meteorology/surface_pressure_met
/Meteorology/temperature_profile_met
/Meteorology/vector_pressure_levels_met
/Meteorology/windspeed_u_met
/Meteorology/windspeed_v_met
/FootprintGeometry/footprint_altitude
/FootprintGeometry/footprint_azimuth
/FootprintGeometry/footprint_land_fraction
/FootprintGeometry/footprint_latitude
/FootprintGeometry/footprint_longitude
/FootprintGeometry/footprint_solar_azimuth
/FootprintGeometry/footprint_solar_zenith
/FootprintGeometry/footprint_stokes_coefficients
/FootprintGeometry/footprint_time
/FootprintGeometry/footprint_time_tai93
/FootprintGeometry/footprint_zenith
/FrameGeometry/relative_velocity
/InstrumentHeader/cnv_coef_highgain_o2
/InstrumentHeader/cnv_coef_highgain_strong_co2
/InstrumentHeader/cnv_coef_highgain_weak_co2
/InstrumentHeader/cnv_coef_medgain_o2
/InstrumentHeader/cnv_coef_medgain_strong_co2
/InstrumentHeader/cnv_coef_medgain_weak_co2
/InstrumentHeader/snr_coef
/Metadata/BuildId
/Metadata/DispersionCoefSamp
/InstrumentHeader/dispersion_coef_samp
/SoundingGeometry
/SoundingGeometry/sounding_id
/SoundingGeometry/sounding_land_water_indicator
/SoundingHeader/gain_swir
/SoundingHeader/sounding_id
/SoundingHeader/wavenumber_coefficients
/SoundingMeasurements
/SoundingMeasurements/radiance_o2
/SoundingMeasurements/radiance_strong_co2
/SoundingMeasurements/radiance_weak_co2
/SoundingSpectra/noise_o2_l1b
/SoundingSpectra/noise_strong_co2_l1b
/SoundingSpectra/noise_weak_co2_l1b
/SoundingSpectra/radiance_o2
/SoundingSpectra/radiance_strong_co2
/SoundingSpectra/radiance_weak_co2
/SpacecraftGeometry/relative_velocity
DATASETS

base_names=""
input_files=`mktemp`
while [ ! -z $1 ] && [ -f $1 ]; do
    echo -n "$1 " >> $input_files
    base_names="${base_names} "`basename $1`
    shift
done
echo >> $input_files

if [ $# -gt 0 ]; then
    sounding_id_file=`mktemp`
    echo $* | sed 's| |\n|g' > ${sounding_id_file}
    # Remove sounding posistion if we have oco files so
    # we are instead using frame ids as splice_product_files.py expects
    if [[ "$base_names" =~ "oco.*\.h5" ]]; then
        sed -i 's|.$||g' $sounding_id_file
    fi

    splice_product_files.py --splice-all -d ${datasets_input_file} -i ${input_files} -s ${sounding_id_file}
    rm -f ${sounding_id_file}
else
    splice_product_files.py --splice-all -d ${datasets_input_file} -i ${input_files} 
fi

rm -f ${input_files}
rm -f ${datasets_input_file}
