begin CASE_SETUP
   remove_existing = True

   begin COPY
      <simulation_runs_dir>/<configuration_case_dir>/<case_spectrum_filename> <case_spectrum_filename>
   end COPY
end CASE_SETUP
