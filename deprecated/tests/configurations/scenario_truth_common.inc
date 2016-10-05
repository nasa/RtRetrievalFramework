# Use soundinginfo for creation of spectra
include = <modifiers_dir>/inp/set_case_soundinginfo.inc
include = <modifiers_dir>/inp/set_spectra_file.inc

# Setup which pressure levels file to use if defined
include = <modifiers_dir>/inp/set_pressure_levels_file.inc

# Set up retrieval vector if user passed in alternative
# Do before aerosol set up since it will modify the retrieval vector and for 
# the scenario truth setup the aerosols are derived from aerosol.dat file
include = <modifiers_dir>/inp/set_retrieval_params.inc

# Use atmosphere, ground and aerosol files from scenario
include = <modifiers_dir>/inp/set_case_atm_apriori.inc
include = <modifiers_dir>/inp/set_case_gnd_apriori.inc
include = <modifiers_dir>/inp/set_case_aer_apriori.inc

# Use static dispersion file
include = <modifiers_dir>/inp/set_static_inst_apriori.inc

# Set up ILS files
include = <modifiers_dir>/inp/set_ils_input.inc

# Set run mode
include = <modifiers_dir>/run/set_run_mode.inc

# Set radiative transfer type
include = <modifiers_dir>/run/set_rt_model.inc

# Turn off or on LSI
include = <modifiers_dir>/run/set_lsi_mode.inc

# Turn on or off polarization if change is defined
include = <modifiers_dir>/run/set_polarization_mode.inc

# Set doppler shift flags
include = <modifiers_dir>/run/set_doppler_switches.inc
