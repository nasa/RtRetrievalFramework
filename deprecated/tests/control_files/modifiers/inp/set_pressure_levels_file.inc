begin FILE
  source_filename = oco_l2.inp

  #####
  # Pressure levels

  begin PICK
      section  = SOUNDING_INFO
      keyword  = pressure_file
      template = <pressure_levels_filename>
      only_if  = 'eval( len("<pressure_levels_filename>") > 0 )'
  end PICK

end FILE


