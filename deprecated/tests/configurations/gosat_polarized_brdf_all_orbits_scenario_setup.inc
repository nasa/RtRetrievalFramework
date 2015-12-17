begin SETTINGS
   # Where to find scenario files
   scenario_directory  = 'eval( os.environ["L2_DATASETS_PATH"] )/scenarios/orbit_simulator/gosat/polarized_brdf'

   # List of scenario IDs from directory listings of orbit directories
   case_list_filename = '<CONTROL_FILE_DIRECTORY>/gosat_all_orbits_land.list'

   # Use orbit?? part of directory when searching for case_ids
   case_dir_levels = 2

   instrument_name = 'gosat'
   brdf_type       = 'lambertian'

end SETTINGS

include = '<CONTROL_FILE_DIRECTORY>/gosat_common_settings.inc'
