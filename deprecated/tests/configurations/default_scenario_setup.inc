begin SETTINGS
   # These should be can overriden

   # List of scenario IDs
   case_list_filename       = '<CONTROL_FILE_DIRECTORY>/orbit_sim_short.list'

   # Where to find scenario files
   scenario_directory       = 'eval( os.environ["L2_DATASETS_PATH"] )/scenarios/orbit_simulator/gosat/pol_brdf_old_aer/orbit_008'

   instrument_name     = 'gosat'
   brdf_type           = 'lambertian'

end SETTINGS

include = '<CONTROL_FILE_DIRECTORY>/gosat_common_settings.inc'
