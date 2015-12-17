begin FILE
  source_filename = oco_l2.inp

  #####
  # All Gases: CO2, O2, H2O

  begin PICK
      section  = PARAMETER_DEFINITION->GAS
      keyword  = a_priori
      template = <case_atmosphere_filename>
  end PICK

  begin PICK
      section  = SOUNDING_INFO
      keyword  = pressure_file
      template = <case_atmosphere_filename>
      not_if   = 'eval( len("<pressure_levels_filename>") > 0 )'
  end PICK

  #############
  # TEMPERATURE

  begin PICK
      section  = PARAMETER_DEFINITION->TEMPERATURE
      keyword  = a_priori
      template = <case_atmosphere_filename>
  end PICK

  #######
  # PSURF

  begin PICK
      section = PARAMETER_DEFINITION->SURFACE_PRESSURE
      keyword = a_priori
      template = <case_psurf_filename>
  end PICK

end FILE


