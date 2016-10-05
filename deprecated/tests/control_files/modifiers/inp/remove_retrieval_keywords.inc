begin FILE
  source_filename = oco_l2.inp

  #######
  # Remove covariance and peturb files with inconsistent sizes
  begin PICK
      section  = PARAMETER_DEFINITION->GAS PARAMETER_DEFINITION->TEMPERATURE PARAMETER_DEFINITION->SURFACE_PRESSURE PARAMETER_DEFINITION->INSTRUMENT PARAMETER_DEFINITION->AEROSOL 
      keyword  = perturb covariance
      delete   = True
  end PICK

  begin PICK
      section  = PARAMETER_DEFINITION->GAS PARAMETER_DEFINITION->TEMPERATURE PARAMETER_DEFINITION->SURFACE_PRESSURE
      keyword  = retrieval_indices scale_parameters
      delete   = True
  end PICK

  begin PICK
      section  = PARAMETER_DEFINITION->BRDF->SPECTRALLY_DEPENDENT
      keyword  = retrieval_indices(1) retrieval_indices(2) retrieval_indices(3)
      delete   = True
      only_if  = 'eval( "<brdf_type>".lower() == "lambertian" )'
  end PICK

  ###########
  # WINDSPEED

  begin PICK
      section  = PARAMETER_DEFINITION->BRDF->SPECTRALLY_INDEPENDENT
      keyword  = retrieval_indices(1) retrieval_indices(2) retrieval_indices(3)
      delete   = True
      only_if  = 'eval( "<brdf_type>".lower() == "coxmunk" )'
  end PICK


  begin PICK
      section  = PARAMETER_DEFINITION->INSTRUMENT
      keyword  = retrieval_indices(1) retrieval_indices(2) retrieval_indices(3)
      delete   = True
      which_section = DISP
  end PICK

end FILE