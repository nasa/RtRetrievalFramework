# Required files
begin CASE_SETUP
   fail_on_error = True
   remove_existing = True

   begin COPY
      <scenario_directory>/<case_id>/atmosphere_<case_short_id>.dat    <case_atmosphere_filename>
      <scenario_directory>/<case_id>/psurf_<case_short_id>.dat         <case_psurf_filename>

      <scenario_directory>/<case_id>/albedo_<case_short_id>.dat        <case_albedo_filename>

      <scenario_directory>/<case_id>/soundinginfo_<case_short_id>.dat  <case_soundinginfo_filename>
   end COPY
end CASE_SETUP

# Optional files
begin CASE_SETUP
   fail_on_error = False
   remove_existing = True
   remove_broken_link = True

   begin LINK
      <scenario_directory>/<case_id>/aerosol_<case_short_id>.dat       <case_aerosol_filename>
      <scenario_directory>/<case_id>/aerosol_od_<case_short_id>.dat    <case_aod_filename>

      <scenario_directory>/<case_id>/windspeed_<case_short_id>.dat     <case_windspeed_filename>
   end LINK
end CASE_SETUP
