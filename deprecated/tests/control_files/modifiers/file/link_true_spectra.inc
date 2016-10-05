begin CASE_SETUP
   link_relative = False
   remove_existing = True

   begin LINK
      <scenario_directory>/<case_id>/spectra_<case_short_id>.dat    <case_spectrum_filename>
   end LINK
end CASE_SETUP
