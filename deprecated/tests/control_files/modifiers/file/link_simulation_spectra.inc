begin CASE_SETUP
   link_relative = True
   fail_on_error = False
   remove_existing = True

   begin LINK
      <simulation_runs_dir>/<configuration_case_dir>/<case_spectrum_filename> <case_spectrum_filename>
   end LINK
end CASE_SETUP
