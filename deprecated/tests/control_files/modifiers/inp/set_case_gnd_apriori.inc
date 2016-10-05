begin FILE
  source_filename = oco_l2.inp

  ########
  # ALBEDO

  begin PICK
      section  = PARAMETER_DEFINITION->BRDF->SPECTRALLY_DEPENDENT
      keyword  = a_priori
      template = <case_albedo_filename>
      only_if  = 'eval( "<brdf_type>".lower() == "lambertian" )'
  end PICK

  ###########
  # WINDSPEED

  begin PICK
      section  = PARAMETER_DEFINITION->BRDF->SPECTRALLY_INDEPENDENT
      keyword  = a_priori
      template = <case_windspeed_filename>
      only_if  = 'eval( "<brdf_type>".lower() == "coxmunk" )'
  end PICK

end FILE


