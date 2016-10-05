# Required files
begin CASE_SETUP
   link_relative = True
   remove_existing = True

   begin MKDIR
      true/
   end MKDIR

   begin LINK
      <simulation_runs_dir>/<configuration_case_dir>/<case_atmosphere_filename>    'true/eval(os.path.basename("<case_atmosphere_filename>"))'
      <simulation_runs_dir>/<configuration_case_dir>/<case_psurf_filename>         'true/eval(os.path.basename("<case_psurf_filename>"))'

      <simulation_runs_dir>/<configuration_case_dir>/<case_albedo_filename>        'true/eval(os.path.basename("<case_albedo_filename>"))'

      <simulation_runs_dir>/<configuration_case_dir>/<case_soundinginfo_filename>  'true/eval(os.path.basename("<case_soundinginfo_filename>"))'
   end LINK
end CASE_SETUP


# Optional files
begin CASE_SETUP
   fail_on_error = False
   remove_existing = True
   remove_broken_link = True

   begin LINK
      <simulation_runs_dir>/<configuration_case_dir>/<case_aerosol_filename>       'true/eval(os.path.basename("<case_aerosol_filename>"))'
      <simulation_runs_dir>/<configuration_case_dir>/<case_aod_filename>           'true/eval(os.path.basename("<case_aod_filename>"))'

      <simulation_runs_dir>/<configuration_case_dir>/<case_windspeed_filename>     'true/eval(os.path.basename("<case_windspeed_filename>"))'
   end LINK
end CASE_SETUP
