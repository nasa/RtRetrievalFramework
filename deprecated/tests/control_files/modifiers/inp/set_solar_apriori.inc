begin FILE
  source_filename = oco_l2.inp

  #######
  # SOLAR

  # Use a different solar line list if so instructed
  begin PICK
      section  = PARAMETER_DEFINITION->SOLAR
      keyword  = solar_linelist
      template = 'eval( os.path.exists("<use_solar_linelist>") and "<use_solar_linelist>" or "./in/static/solar/<use_solar_linelist>" )'
      only_if  = 'eval( len("<use_solar_linelist>") > 0 )'
  end PICK

  # Use a different solar continuum parameterization if so instructed
  begin PICK
      section  = PARAMETER_DEFINITION->SOLAR
      keyword  = a_priori
      template = 'eval( os.path.exists("<use_solar_cont>") and "<use_solar_cont>" or "./in/static/solar/<use_solar_cont>" )'
      only_if  = 'eval( len("<use_solar_cont>") > 0 )'
  end PICK


end FILE
