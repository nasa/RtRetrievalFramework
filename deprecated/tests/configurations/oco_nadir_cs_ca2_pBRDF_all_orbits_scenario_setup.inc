begin SETTINGS
   # Where to find scenario files
   scenario_directory  = 'eval( os.environ["L2_DATASETS_PATH"] )/scenarios/orbit_simulator/oco/nadir_cs_ca2_pBRDF'

   # List of scenario IDs from directory listings of orbit directories
   case_list_filename = '<CONTROL_FILE_DIRECTORY>/oco_nadir_cs_ca2_pBRDF_all_orbits.list'

   # Use orbit?? part of directory when searching for case_ids
   case_dir_levels = 2

   instrument_name = 'oco'
   brdf_type       = 'lambertian'

end SETTINGS

include = '<CONTROL_FILE_DIRECTORY>/oco_common_settings.inc'
