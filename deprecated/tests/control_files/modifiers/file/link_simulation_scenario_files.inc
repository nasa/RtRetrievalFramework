# Required files
begin CASE_SETUP
   fail_on_error = True
   remove_existing = True
   link_relative = True

   begin LINK
      <simulation_runs_dir>/<configuration_case_dir>/<case_atmosphere_filename>    <case_atmosphere_filename>
      <simulation_runs_dir>/<configuration_case_dir>/<case_psurf_filename>         <case_psurf_filename>

      <simulation_runs_dir>/<configuration_case_dir>/<case_albedo_filename>        <case_albedo_filename>

      <simulation_runs_dir>/<configuration_case_dir>/<case_soundinginfo_filename>  <case_soundinginfo_filename>
   end LINK
end CASE_SETUP


# Optional files
begin CASE_SETUP
   fail_on_error = True
   remove_existing = True
   remove_broken_link = True

   begin LINK
      <simulation_runs_dir>/<configuration_case_dir>/<case_aerosol_filename>       <case_aerosol_filename>
      <simulation_runs_dir>/<configuration_case_dir>/<case_aod_filename>           <case_aod_filename>

      <simulation_runs_dir>/<configuration_case_dir>/<case_windspeed_filename>     <case_windspeed_filename>
   end LINK
end CASE_SETUP
