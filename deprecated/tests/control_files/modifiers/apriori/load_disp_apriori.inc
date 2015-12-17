begin CASE_SETUP
   remove_existing = True
   only_if = 'eval( "<instrument_name>".lower() == "gosat" )'

   begin MKDIR
      in/instrument
   end MKDIR
end CASE_SETUP

begin FILE
   source_filename = <case_spectrum_filename>
   dest_filename   = <case_dispersion_filename>

   begin LOAD_DISP_APRIORI
      dispersion_in_file = in/static/instrument/dispersion.dat
      only_if = 'eval( "<instrument_name>".lower() == "gosat" )'
      index_scheme = '<disp_index_scheme|>'
   end LOAD_DISP_APRIORI
end FILE

begin FILE
   source_filename = oco_l2.inp

   begin PICK 
      section  = PARAMETER_DEFINITION->INSTRUMENT
      keyword  = a_priori
      which_section = DISP
      template = <case_dispersion_filename>
      only_if = 'eval( "<instrument_name>".lower() == "gosat" )'
   end PICK 

end FILE
