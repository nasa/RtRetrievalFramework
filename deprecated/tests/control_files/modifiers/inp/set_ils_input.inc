begin KEYWORDS
   # Always set default for oco instrument but not for gosat
   ils_type_name = 'eval( "<instrument_name>".lower() == "oco" and "<ils_type_name|table>" or "<ils_type_name|>" )'
end KEYWORDS

begin MAP
   name = ils_sinc_inst_map
   begin VALUES
      oco    ''
      gosat  in/static/instrument/ils_sinc.dat
   end VALUES
end MAP

begin MAP
   name = ils_gauss_inst_map
   begin VALUES
      oco    in/static/instrument/ils_gauss.dat
      gosat  ''
   end VALUES
end MAP

begin MAP
   name = ils_table_inst_map
   begin VALUES
      oco        <l2_input_dir>/oco/ils_tables/ils_table_orbit_sim.dat
      oco_tconn  <l2_input_dir>/oco/ils_tables/ils_table_5.dat
      gosat      in/static/instrument/TFTSILSF_A006P006_B1S.dat in/static/instrument/TFTSILSF_A006P006_B2.dat in/static/instrument/TFTSILSF_A006P006_B3.dat
   end VALUES
end MAP

begin MAP
   name = ils_type_map
   begin VALUES
      sinc   '[ils_sinc_inst_map:instrument_name]'
      gauss  '[ils_gauss_inst_map:instrument_name]'
      table  '[ils_table_inst_map:instrument_name]'
   end VALUES
end MAP

not_if = 'eval(len("<ils_type_name>") > 0 and len("[ils_type_map:ils_type_name]") <= 0)'

begin FILE
  source_filename = oco_l2.inp

  ######
  # ILS
  begin PICK
      only_if = 'eval(len("<ils_type_name>") > 0)'
      section = PARAMETER_DEFINITION->INSTRUMENT
      keyword = 'a_priori(1)' 'a_priori(2)' 'a_priori(3)'
      delete = True
      which_section = ILS
  end PICK

  begin PICK
      only_if = 'eval(len("<ils_type_name>") > 0)'
      section = PARAMETER_DEFINITION->INSTRUMENT
      keyword = a_priori
      template = [ils_type_map:ils_type_name]
      which_section = ILS
      index_format = '{keyword}({index})'
  end PICK

end FILE

