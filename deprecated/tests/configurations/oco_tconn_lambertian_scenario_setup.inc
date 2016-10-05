begin SETTINGS
   # Where to find scenario files
   scenario_directory  = 'eval( os.environ["L2_DATASETS_PATH"] )/scenarios/oco_tconn'

   # List of scenario IDs from directory listings of orbit directories
   case_list_filename = '<CONTROL_FILE_DIRECTORY>/oco_tconn_lambertian_cases.list'

   # Use orbit?? part of directory when searching for case_ids
   case_dir_levels = 2

   instrument_name = 'oco_tconn'
   brdf_type       = 'lambertian'

   ils_type_name   = 'table'
   use_original_aerosol_mode = True
end SETTINGS

include = '<CONTROL_FILE_DIRECTORY>/oco_tconn_common_settings.inc'
